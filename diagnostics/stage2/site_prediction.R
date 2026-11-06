# Patient-grouped local cross-validation, not external-site outcome validation.
source('diagnostics/stage2/lib.R');source('diagnostics/pooled_roc/lib.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=3L)stop('Usage: site_prediction.R SITE BUNDLE NEW_OUTPUT')
site<-a[1];bundle<-a[2];out<-a[3];if(dir.exists(out))stop('Choose fresh prediction output')
plan<-jsonlite::fromJSON(file.path(bundle,'analysis_plan.json'))
d<-ap_data(site,'full');d<-d[d$valid_outcome,,drop=FALSE]
# In-hospital endpoints have complete ascertainment through exit or day 28.
stopifnot(all(d$time_from_landmark>0))
patients<-sort(unique(d$patient_id),method='radix');B<-plan$prediction$bootstrap
set.seed(plan$prediction$seed+sum(utf8ToInt(site)));P<-length(patients)
pw<-if(P)replicate(B,tabulate(sample.int(P,P,replace=TRUE),nbins=P)) else matrix(0,0,B)
cohorts<-c('all_case_missing_indicators','measured_SOFA','complete_hb_proxy','matched_complete_SOFA')
rows<-list();metrics<-list();den<-list();foldchecks<-list();coeff<-list();calibration<-list();decisions<-list();private<-list()
for(cohort in cohorts){
 cat(site,'local CV cohort:',cohort,'\n');flush.console()
 z<-if(cohort=='measured_SOFA')d[!is.na(d$sofa_total),] else ap_cases(d,cohort)
 kinds<-if(cohort=='measured_SOFA')c('baseline','phenotype','sofa_only','baseline_sofa','summaries','summaries_phenotype','sofa_summaries','sofa_summaries_phenotype') else ap_models(cohort)
 fold<-s2_folds(z$patient_id,plan$prediction$folds,plan$prediction$seed)
 for(endpoint in plan$prediction$endpoints){
  y<-as.integer(z$time_from_landmark<=28 & if(endpoint=='death')z$event_type=='death' else z$event_type%in%c('death','hospice'))
  preds<-list();audits<-list()
  for(kind in kinds){
   if(!nrow(z)){preds[[kind]]<-numeric();next}
   x<-ap_matrix(z,kind,if(cohort=='measured_SOFA')'complete_hb_proxy' else cohort,ap_phenotypes('full'))
   f<-s2_oof(x,y,z$patient_id,fold);preds[[kind]]<-f$prediction
   audits[[kind]]<-f$audit |> dplyr::mutate(site=site,cohort=cohort,endpoint=endpoint,model=kind)
   fit<-if(length(unique(y))==2L)tryCatch(ridge_logistic(x,y,z$patient_id),error=identity) else NULL
   if(is.list(fit)&&!inherits(fit,'error')&&isTRUE(fit$converged))coeff[[paste(cohort,endpoint,kind)]]<-tibble::tibble(site=site,cohort=cohort,endpoint=endpoint,model=kind,term=names(fit$beta),estimate=unname(fit$beta),std_error=sqrt(pmax(0,diag(fit$covariance))),status='ok')
  }
  foldchecks[[paste(cohort,endpoint)]]<-dplyr::bind_rows(audits)
  # Shared evaluation patients across every model in this cohort.
  usable<-if(nrow(z))Reduce(`&`,lapply(preds,is.finite)) else logical()
  q<-z[usable,,drop=FALSE];weights<-pw[match(q$patient_id,patients),,drop=FALSE]
  den[[paste(cohort,endpoint)]]<-tibble::tibble(site=site,cohort=cohort,endpoint=endpoint,candidate_n=nrow(z),n=nrow(q),patients=dplyr::n_distinct(q$patient_id),excluded_nonestimable_fold_n=sum(!usable),status=if(nrow(q))'ok' else 'nonestimable_or_no_complete_cases')
  for(kind in kinds){
   p<-preds[[kind]][usable];bins<-roc_bin(p)
   target<-if(endpoint=='death')q$event_type=='death' else q$event_type%in%c('death','hospice')
   case_day<-ifelse(target & q$time_from_landmark<=28,pmax(1,ceiling(q$time_from_landmark)),29)
   groups<-paste(bins,case_day,sep='|');keys<-unique(groups)
   sh<-tibble::tibble(bin=as.integer(sub('[|].*','',keys)),case_day=as.integer(sub('.*[|]','',keys)),n=tabulate(match(groups,keys),nbins=length(keys)))
   for(b in seq_len(B)){
    ns<-numeric(length(keys));if(length(keys)){rr<-rowsum(weights[,b],match(groups,keys),reorder=FALSE);ns[as.integer(rownames(rr))]<-rr[,1]};sh[[sprintf('n_b%03d',b)]]<-ns
   }
   if(!nrow(sh)){sh<-tibble::tibble(bin=0L,case_day=29L,n=0);for(b in seq_len(B))sh[[sprintf('n_b%03d',b)]]<-0}
   rows[[paste(cohort,endpoint,kind)]]<-sh |> dplyr::mutate(site=site,cohort=cohort,endpoint=endpoint,model=kind)
   for(day in plan$prediction$tdauc_days){
    event<-as.integer(q$time_from_landmark<=day & if(endpoint=='death')q$event_type=='death' else q$event_type%in%c('death','hospice'))
    h<-roc_hist(bins,event)
    # Always emit an empty bin for nonestimable cohorts, with zero denominators.
    if(!length(h$bin))h<-list(bin=0L,cases=0,controls=0)
    bs<-matrix(NA_real_,B,2)
    for(b in seq_len(B)){
     hb<-roc_hist(bins,event,weights[,b],h$bin)
     bs[b,]<-c(roc_auc(hb$cases,hb$controls),if(length(event)&&sum(weights[,b])>0)weighted.mean((event-p)^2,weights[,b]) else NA_real_)
    }
    exact<-auc_rank(event,p);approx<-roc_auc(h$cases,h$controls)
    bound<-if(sum(event)&&sum(1-event))sum(h$cases*h$controls)/(2*sum(h$cases)*sum(h$controls)) else NA_real_
    if(is.finite(exact))stopifnot(abs(exact-approx)<=bound+1e-12)
    ci<-if(any(is.finite(bs[,1])))quantile(bs[,1],c(.025,.975),na.rm=TRUE) else c(NA,NA)
    metrics[[paste(cohort,endpoint,kind,day)]]<-tibble::tibble(site=site,cohort=cohort,endpoint=endpoint,model=kind,day=day,n=nrow(q),cases=sum(event),controls=sum(1-event),exact_local_auc=exact,binned_local_auc=approx,absolute_difference=abs(exact-approx),auc_quantization_bound=bound,conf_low=ci[1],conf_high=ci[2],brier=if(length(event))mean((event-p)^2) else NA_real_,status=if(is.finite(exact))'ok' else 'nonestimable_endpoint_partition')
   }
   event<-as.integer(q$time_from_landmark<=28 & if(endpoint=='death')q$event_type=='death' else q$event_type%in%c('death','hospice'))
   if(length(p)){
    cb<-tibble::tibble(bin=pmin(9L,floor(p*10)),p=p,event=event) |> dplyr::group_by(bin) |> dplyr::summarise(n=dplyr::n(),predicted_sum=sum(p),events=sum(event),.groups='drop') |> dplyr::mutate(site=site,cohort=cohort,endpoint=endpoint,model=kind)
    calibration[[paste(cohort,endpoint,kind)]]<-cb
    for(t in plan$prediction$decision_thresholds)decisions[[paste(cohort,endpoint,kind,t)]]<-tibble::tibble(site=site,cohort=cohort,endpoint=endpoint,model=kind,threshold=t,n=length(event),true_positive=sum(event*(p>=t)),false_positive=sum((1-event)*(p>=t)),net_benefit=mean(event*(p>=t))-mean((1-event)*(p>=t))*t/(1-t))
   }
  }
  private[[paste(cohort,endpoint)]]<-list(ids=z$hospitalization_id,patient_id=z$patient_id,fold=fold,prediction=preds,usable=usable)
 }
}
dir.create(out,recursive=TRUE)
for(pair in list(list(x=rows,name='score_histograms'),list(x=metrics,name='local_validation'),list(x=den,name='cohort_counts'),list(x=foldchecks,name='fold_audit'),list(x=coeff,name='coefficients'),list(x=calibration,name='calibration_bins'),list(x=decisions,name='decision_curves'))){x<-dplyr::bind_rows(pair$x);if(!ncol(x))x<-tibble::tibble(site=character(),status=character());fed_write(x,file.path(out,paste0(pair$name,'.csv')))}
dir.create(file.path(ap_root(),'sites',site,'local'),recursive=TRUE,showWarnings=FALSE);saveRDS(private,file.path(ap_root(),'sites',site,'local/internal_validation.rds'))
files<-list.files(out,full.names=TRUE)
fed_json(list(status='passed',site=site,consensus_id=jsonlite::fromJSON(file.path(ap_primary(),'coordinator/full/frozen_manifest.json'))$consensus_id,validation='patient_grouped_local_cross_validation',external_site_outcome_validation=FALSE,folds=plan$prediction$folds,bins=ROC_BINS,bootstrap=B,bootstrap_scope='patient-cluster resampling of out-of-fold predictions; models not refit in bootstrap; conditional on frozen phenotypes and fitted folds',histogram='Counts by fixed probability bin and first daily event horizon; case_day=29 denotes no target event through day28',score='Out-of-fold 28-day hospital endpoint risk evaluated over days 1-28',control_definition='Cumulative/dynamic AUC2: competing exits are controls; no ordinary censoring before 28 days among valid-outcome hospitalizations',cohorts=cohorts,models='See cohort/model rows; four-model comparisons use matched_complete_SOFA',unknown_outcomes_excluded_from_prediction=sum(!ap_data(site,'full')$valid_outcome),files=as.list(setNames(unname(tools::md5sum(files)),basename(files)))),file.path(out,'manifest.json'))
