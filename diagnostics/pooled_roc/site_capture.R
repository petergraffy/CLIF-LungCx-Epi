source('diagnostics/apache_adapted/lib.R');source('diagnostics/pooled_roc/lib.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=4)stop('Usage: site_capture.R SITE PRIMARY APACHE NEW_SHARE')
site<-a[1];primary<-a[2];apache<-a[3];out<-a[4];if(dir.exists(out))stop('Choose fresh output')
manifest<-jsonlite::fromJSON(file.path(primary,'coordinator/full/frozen_manifest.json'))
full<-readRDS(file.path(primary,'sites',site,'local',paste0('loo_',site),'outcomes.rds'))$data
cache<-readRDS(file.path(apache,'sites',site,'local/evaluation.rds'))
transport<-jsonlite::fromJSON(file.path(apache,'sites',site,'export/transport_manifest.json'))
stopifnot(transport$site==site,!site%in%transport$training_sites,!transport$held_out_refitting)
models<-c('baseline','baseline_sofa','baseline_adapted','phenotype');B<-200L
rows<-list();checks<-list();cohorts<-list();ids_reference<-NULL;weights<-NULL
for(endpoint in c('death_or_hospice','death')){
 x<-cache[[paste('matched_complete_SOFA',endpoint)]]
 d<-full[match(x$ids,full$hospitalization_id),]
 stopifnot(identical(d$hospitalization_id,x$ids),identical(d$patient_id,x$patient_id),all(d$valid_outcome),all(d$time_from_landmark>0))
 if(is.null(ids_reference)){
  ids_reference<-x$ids;patients<-unique(d$patient_id);P<-length(patients);if(P==0)stop('No matched complete-score evaluation patients; four-model comparison unavailable')
  set.seed(20261101+sum(utf8ToInt(site)));pw<-replicate(B,tabulate(sample.int(P,P,replace=TRUE),nbins=P))
  weights<-pw[match(d$patient_id,patients),,drop=FALSE]
 }else stopifnot(identical(ids_reference,x$ids))
 for(model in models){
  p<-x$pred[[model]];stopifnot(length(p)==nrow(d),all(is.finite(p)))
  bins<-roc_bin(p)
  for(day in c(7,14,28)){
   event<-as.integer(d$time_from_landmark<=day & if(endpoint=='death')d$event_type=='death' else d$event_type%in%c('death','hospice'))
   h<-roc_hist(bins,event)
   z<-tibble::tibble(site=site,endpoint=endpoint,model=model,day=day,bin=h$bin,cases=h$cases,controls=h$controls)
   for(b in seq_len(B)){
    q<-roc_hist(bins,event,weights[,b],h$bin)
    z[[sprintf('cases_b%03d',b)]]<-q$cases;z[[sprintf('controls_b%03d',b)]]<-q$controls
   }
   rows[[paste(endpoint,model,day)]]<-z
   exact<-auc_rank(event,p);approx<-roc_auc(h$cases,h$controls)
   bound<-if(sum(event)&&sum(1-event))sum(h$cases*h$controls)/(2*sum(h$cases)*sum(h$controls)) else NA_real_
   if(is.finite(exact))stopifnot(abs(exact-approx)<=bound+1e-12)
   checks[[paste(endpoint,model,day)]]<-tibble::tibble(site=site,endpoint=endpoint,model=model,day=day,n=nrow(d),cases=sum(event),controls=sum(1-event),exact_local_auc=exact,binned_local_auc=approx,absolute_difference=abs(exact-approx),auc_quantization_bound=bound)
  }
 }
 cohorts[[endpoint]]<-tibble::tibble(site=site,endpoint=endpoint,n=nrow(d),patients=length(unique(d$patient_id)),cohort='complete modified SOFA and complete adapted APACHE; identical cases across all four models')
}
dir.create(out,recursive=TRUE)
for(pair in list(list(x=dplyr::bind_rows(rows),name='score_histograms.csv'),list(x=dplyr::bind_rows(checks),name='local_validation.csv'),list(x=dplyr::bind_rows(cohorts),name='cohort_counts.csv'))){
 stopifnot(!any(names(pair$x)%in%c('patient_id','hospitalization_id','event_time','t0','marker','prediction')))
 readr::write_csv(pair$x,file.path(out,pair$name),na='')
}
f<-list.files(out,full.names=FALSE)
jsonlite::write_json(list(status='passed',site=site,consensus_id=manifest$consensus_id,bins=ROC_BINS,bootstrap=B,models=models,days=c(7,14,28),endpoints=c('death','death_or_hospice'),control_definition='timeROC AUC2: all noncases including competing exits; no ordinary censoring in evaluation cohort',score='Frozen held-out 28-day model risk, evaluated at each horizon; common probability bins of width 0.00001',bootstrap_unit='Within-site patients; all admissions retained; shared weights across all models/horizons/endpoints',held_out_refitting=FALSE,training_sites=transport$training_sites,files=as.list(setNames(unname(tools::md5sum(file.path(out,f))),f))),file.path(out,'manifest.json'),pretty=TRUE,auto_unbox=TRUE)
cat(site,'four-model matched-cohort ROC aggregate histograms exported\n')
