source('diagnostics/apache_adapted/lib.R');site<-commandArgs(trailingOnly=TRUE)[1];version<-paste0('loo_',site)
home<-file.path(ap_root(),'sites',site);out<-file.path(home,'export');d<-ap_data(site,version)
predpath<-file.path(ap_root(),'coordinator',version,'coefficients.csv');b<-readr::read_csv(predpath,show_col_types=FALSE);m<-jsonlite::fromJSON(file.path(dirname(predpath),'frozen_predictor.json'))
stopifnot(!site%in%m$training_sites,m$spec_hash==fed_hash(ap_spec()),m$predictor_hash==digest::digest(file=predpath,algo='sha256'))
rows<-list();paired<-list();decision<-list();bins<-list();private<-list();metadata<-list()
q95<-function(x)if(sum(is.finite(x))>=100)unname(quantile(x[is.finite(x)],c(.025,.975))) else c(NA_real_,NA_real_)
for(cohort in ap_spec()$cohorts){z<-ap_cases(d,cohort)
 for(endpoint in ap_spec()$endpoints){
 y<-as.integer(z$time_from_landmark<=28 & if(endpoint=='death')z$event_type=='death' else z$event_type%in%c('death','hospice'))
 safe<-sum(y)>=aggregate_min_cell(10)&&sum(1-y)>=aggregate_min_cell(10);pred<-list();metrics<-list();status<-list();thresholds<-seq(.05,.6,.05)
 for(kind in ap_models(cohort)){
  bb<-b[b$cohort==cohort&b$endpoint==endpoint&b$model==kind,]
  required<-c('(Intercept)',if(grepl('adapted',kind))'adapted_APS_age',if(grepl('sofa',kind))'modified_SOFA')
  transportable<-all(required%in%bb$term)&&all(bb$contributing_sites[match(required,bb$term)]>0)
  if(!nrow(z)||!transportable){pred[[kind]]<-rep(NA_real_,nrow(z));metrics[[kind]]<-setNames(rep(NA_real_,4),c('auc','brier','calibration_intercept','calibration_slope'));status[[kind]]<-if(!nrow(z))'no_evaluation_cases' else 'required_predictive_terms_unestimated'} else {
   x<-ap_matrix(z,kind,cohort,ap_phenotypes(version));stopifnot(setequal(colnames(x),bb$term));beta<-setNames(bb$estimate,bb$term)
   pred[[kind]]<-plogis(drop(x%*%beta[colnames(x)]));metrics[[kind]]<-pred_metrics(y,pred[[kind]]);status[[kind]]<-if(safe)'ok' else'withheld_small_endpoint_partition'
  }
 }
 pairs<-ap_pairs(cohort);boot<-list();dc<-list();by<-split(seq_len(nrow(z)),z$patient_id);patients<-names(by)
 set.seed(20261030)
 if(length(patients))for(i in 1:200){
  ix<-unlist(by[sample(patients,length(patients),replace=TRUE)],use.names=FALSE)
  bm<-lapply(pred,function(pr)if(anyNA(pr))c(auc=NA_real_,brier=NA_real_) else c(auc=auc_rank(y[ix],pr[ix]),brier=mean((y[ix]-pr[ix])^2)))
  vals<-unlist(lapply(names(bm),function(k)setNames(bm[[k]],paste0(k,'_',names(bm[[k]])))))
  for(pair in names(pairs))for(metric in c('auc','brier'))vals[paste(pair,metric,sep='_')]<-bm[[pairs[[pair]][1]]][metric]-bm[[pairs[[pair]][2]]][metric]
  boot[[i]]<-vals
  dc[[i]]<-unlist(lapply(pred,function(pr)vapply(thresholds,function(t)mean(y[ix]*(pr[ix]>=t))-mean((1-y[ix])*(pr[ix]>=t))*t/(1-t),numeric(1))))
 }
 boots<-if(length(boot))do.call(rbind,boot) else NULL;db<-if(length(dc))do.call(rbind,dc) else NULL
 for(kind in names(metrics)){
  ok<-safe&&status[[kind]]=='ok'
  for(metric in names(metrics[[kind]])){
   ci<-if(!is.null(boots)&&metric%in%c('auc','brier'))q95(boots[,paste0(kind,'_',metric)]) else c(NA_real_,NA_real_)
   rows[[paste(cohort,endpoint,kind,metric)]]<-tibble::tibble(site=site,cohort=cohort,endpoint=endpoint,model=kind,metric=metric,estimate=if(ok)unname(metrics[[kind]][metric]) else NA_real_,conf_low=if(ok)ci[1] else NA_real_,conf_high=if(ok)ci[2] else NA_real_,status=status[[kind]])
  }
  for(j in seq_along(thresholds)){
   t<-thresholds[j];v<-pred[[kind]]>=t;ci<-if(!is.null(db))q95(db[,(match(kind,names(pred))-1)*length(thresholds)+j]) else c(NA_real_,NA_real_)
   decision[[paste(cohort,endpoint,kind,j)]]<-tibble::tibble(site=site,cohort=cohort,endpoint=endpoint,model=kind,threshold=t,net_benefit=if(ok)mean(y*v)-mean((1-y)*v)*t/(1-t) else NA_real_,conf_low=if(ok)ci[1] else NA_real_,conf_high=if(ok)ci[2] else NA_real_,treat_all=if(ok)mean(y)-(1-mean(y))*t/(1-t) else NA_real_,treat_none=0,status=status[[kind]])
  }
  if(nrow(z)&&all(is.finite(pred[[kind]]))){
   cal<-data.frame(bin=cut(pred[[kind]],breaks=seq(0,1,.2),include.lowest=TRUE),y=y,p=pred[[kind]]) |> dplyr::group_by(bin) |> dplyr::summarise(n=dplyr::n(),events_n=sum(y),non_events_n=sum(1-y),mean_predicted=mean(p),observed_risk=mean(y),.groups='drop')
   cal<-protect_aggregate(cal,c('n','events_n','non_events_n'),character(),10);cal$site<-site;cal$cohort<-cohort;cal$endpoint<-endpoint;cal$model<-kind
   bins[[paste(cohort,endpoint,kind)]]<-cal
  }
 }
 for(pair in names(pairs))for(metric in c('auc','brier')){
  ok<-safe&&all(vapply(pairs[[pair]],function(k)status[[k]]=='ok',logical(1)))
  ci<-if(!is.null(boots))q95(boots[,paste(pair,metric,sep='_')]) else c(NA_real_,NA_real_)
  delta<-metrics[[pairs[[pair]][1]]][metric]-metrics[[pairs[[pair]][2]]][metric]
  paired[[paste(cohort,endpoint,pair,metric)]]<-tibble::tibble(site=site,cohort=cohort,endpoint=endpoint,comparison=pair,metric=metric,estimate=if(ok)unname(delta) else NA_real_,conf_low=if(ok)ci[1] else NA_real_,conf_high=if(ok)ci[2] else NA_real_,status=if(ok)'ok' else'unavailable_or_withheld')
 }
 nb<-fed_band_n(nrow(z));eb<-fed_band_n(sum(y));metadata[[paste(cohort,endpoint)]]<-tibble::tibble(site=site,cohort=cohort,endpoint=endpoint,n_lower=nb[1],n_upper=nb[2],events_lower=eb[1],events_upper=eb[2],same_cases_across_models=TRUE,status=if(safe)'ok' else'withheld_small_endpoint_partition')
 private[[paste(cohort,endpoint)]]<-list(ids=z$hospitalization_id,patient_id=z$patient_id,y=y,pred=pred,bootstrap=boots,n=nrow(z),events=sum(y))
 }
}
for(pair in list(list(x=rows,file='transport_metrics.csv'),list(x=paired,file='paired_increments.csv'),list(x=decision,file='decision_curves.csv'),list(x=bins,file='calibration_bins.csv'),list(x=metadata,file='evaluation_metadata.csv')))ap_export(dplyr::bind_rows(pair$x),file.path(out,pair$file))
saveRDS(private,file.path(home,'local/evaluation.rds'))
fed_json(list(site=site,held_out_site=site,training_sites=m$training_sites,held_out_refitting=FALSE,bootstrap=200,bootstrap_scope='patient bootstrap conditional on frozen learner; calibration point estimates',same_cases_within_cohort=TRUE,specification=ap_spec()),file.path(out,'transport_manifest.json'))
cat(site,'held-out adapted-score evaluation complete\n')
