# Local-only survival and time-dependent discrimination calculations.
km_release<-function(d,time,group='Overall',days=0:28,min_cell=10L){
 min_cell<-aggregate_min_cell(min_cell)
 out<-list()
 for(endpoint in c('death','death_or_hospice')){
  event<-if(endpoint=='death')d$valid_outcome & d$event_type=='death' else d$valid_outcome & d$event_type%in%c('death','hospice')
  fit<-survival::survfit(survival::Surv(time,event)~1,id=d$hospitalization_id,cluster=d$patient_id,robust=TRUE,conf.type='log-log')
  grid<-days[days<=max(time)];sm<-summary(fit,times=grid)
  nb<-aggregate_count_band(nrow(d));rb<-t(vapply(sm$n.risk,aggregate_count_band,numeric(2)))
  z<-tibble::tibble(group=group,endpoint=endpoint,day=sm$time,survival=sm$surv,std_error=sm$std.err,conf_low=sm$lower,conf_high=sm$upper,
                    n=nrow(d),n_at_risk=sm$n.risk,n_lower=nb[1],n_upper=nb[2],
                    at_risk_lower=rb[,1],at_risk_upper=rb[,2],cell_status='observed')
  for(j in seq_len(nrow(z))){t<-z$day[j];counts<-c(sum(event & time<=t),sum(!event & time<=t),sm$n.risk[j])
   if(any(counts>0 & counts<min_cell)){z[j,c('survival','std_error','conf_low','conf_high')]<-NA_real_;z$cell_status[j]<-'suppressed_partition'}
  }
  out[[endpoint]]<-z
 }
 dplyr::bind_rows(out)
}

td_auc_rank<-function(T,delta,marker,t,control=2L){
 cases<-T<=t & delta==1L
 controls<-if(control==1L)T>t else !cases
 # This exact simplification is valid only with no ordinary right censoring.
 stopifnot(all(delta%in%c(1L,2L)))
 y<-as.integer(cases[cases|controls]);m<-marker[cases|controls]
 n1<-sum(y);n0<-length(y)-n1
 if(n1==0||n0==0)return(NA_real_)
 (sum(rank(m,ties.method='average')[y==1])-n1*(n1+1)/2)/(n1*n0)
}

td_capture<-function(d,predictions,cohort,endpoint,site,bootstrap=200L,min_cell=10L){
 min_cell<-max(1L,aggregate_min_cell(min_cell))
 T<-d$time_from_landmark
 delta<-if(endpoint=='death')ifelse(d$event_type=='death',1L,2L) else ifelse(d$event_type%in%c('death','hospice'),1L,2L)
 stopifnot(all(d$valid_outcome),all(T>0),all(d$event_type%in%c('death','hospice','live_discharge')))
 stopifnot(all(vapply(predictions,function(p)length(p)==nrow(d)&&all(is.finite(p)),logical(1))))
 days<-1:28;models<-names(predictions);points<-list();rocs<-list();reference<-list()
 for(model in models){
  # Tiny offset implements inclusive <= daily horizons despite timeROC's strict <.
  q<-timeROC::timeROC(T=T,delta=delta,marker=predictions[[model]],cause=1,times=days+1e-8,weighting='marginal',ROC=TRUE,iid=FALSE)
  w<-q$weights$IPCW.times
  # Beyond the last hospital exit, event-free controls are empty. timeROC sets
  # G(t) to NA there, while AUC_2 remains estimable from known competing exits.
  stopifnot(all(w[is.finite(w)]==1),all(q$times[!is.finite(w)]>=max(T)))
  reference[[model]]<-q
  for(control in 1:2){
   values<-vapply(days,function(t)td_auc_rank(T,delta,predictions[[model]],t,control),numeric(1))
   official<-if(control==1)q$AUC_1 else q$AUC_2
   stopifnot(all(abs(values-official)<1e-10,na.rm=TRUE))
  }
  for(t in c(7,14,28))for(control in 1:2){
   x<-if(control==1)q$FP_1[,t] else q$FP_2[,t];y<-q$TP[,t]
   good<-is.finite(x)&is.finite(y)
   if(sum(good)>1){
    # Fixed FPR grid returns no score cutoffs or patient-specific ROC vertices.
    xx<-x[good];yy<-y[good];u<-sort(unique(xx));top<-vapply(u,function(v)max(yy[xx==v]),numeric(1))
    fpr<-seq(0,1,.02);sens<-approx(u,top,xout=fpr,rule=2,ties='ordered')$y
    cases<-sum(T<=t&delta==1);controls<-if(control==1)sum(T>t) else nrow(d)-cases
    safe<-cases>=min_cell&&controls>=min_cell
    rocs[[paste(model,t,control)]]<-tibble::tibble(site=site,cohort=cohort,endpoint=endpoint,model=model,day=t,control_definition=control,fpr=fpr,sensitivity=if(safe)sens else NA_real_,auc=if(safe)(if(control==1)q$AUC_1[t] else q$AUC_2[t]) else NA_real_,status=if(safe)'observed' else if(aggregate_exact())'no_cases_or_controls' else 'withheld_small_case_control_partition')
   }
  }
 }
 # Same patient resamples across models, days and control definitions.
 by<-split(seq_len(nrow(d)),d$patient_id);set.seed(20261030)
 boots<-array(NA_real_,dim=c(bootstrap,length(models),28L,2L),dimnames=list(NULL,models,NULL,NULL))
 for(b in seq_len(bootstrap)){
  ix<-unlist(by[sample(names(by),length(by),replace=TRUE)],use.names=FALSE)
  for(j in seq_along(models))for(control in 1:2)boots[b,j,,control]<-vapply(days,function(t)td_auc_rank(T[ix],delta[ix],predictions[[models[j]]][ix],t,control),numeric(1))
 }
 for(j in seq_along(models))for(t in days)for(control in 1:2){
  model<-models[j];case_n<-sum(T<=t&delta==1);control_n<-if(control==1)sum(T>t) else nrow(d)-case_n
  safe<-case_n>=min_cell&&control_n>=min_cell;v<-boots[,j,t,control];ci<-if(sum(is.finite(v))>=100)quantile(v,c(.025,.975),na.rm=TRUE) else c(NA_real_,NA_real_)
  est<-if(control==1)reference[[model]]$AUC_1[t] else reference[[model]]$AUC_2[t]
  cb<-aggregate_count_band(case_n);nb<-aggregate_count_band(control_n)
  points[[paste(model,t,control)]]<-tibble::tibble(site=site,cohort=cohort,endpoint=endpoint,model=model,day=t,control_definition=control,estimate=if(safe)est else NA_real_,conf_low=if(safe)ci[1] else NA_real_,conf_high=if(safe)ci[2] else NA_real_,cases=case_n,controls=control_n,case_lower=cb[1],case_upper=cb[2],control_lower=nb[1],control_upper=nb[2],bootstrap_replicates=sum(is.finite(v)),status=if(safe)'observed' else if(aggregate_exact())'no_cases_or_controls' else 'withheld_small_case_control_partition')
 }
 list(metrics=dplyr::bind_rows(points),roc=dplyr::bind_rows(rocs))
}
