source('diagnostics/revision20261022/score_lib.R')
site<-commandArgs(trailingOnly=TRUE)[1];home<-file.path(rev_root(),'sites',site);out<-file.path(home,'export','scores');rows<-list();meta<-list();fits<-list()
for(version in paste0('loo_',c('UCMC','NU','RUSH'))){
 d<-readRDS(file.path(old_root(),'sites',site,'local',version,'outcomes.rds'))$data
 for(cohort in score_spec()$cohorts){
 z<-if(cohort=='measured_sofa')d[!is.na(d$sofa_total),] else d
 for(endpoint in score_spec()$endpoints){
 y<-as.integer(z$time_from_landmark<=28 & if(endpoint=='death')z$event_type=='death' else z$event_type%in%c('death','hospice'))
 for(kind in score_spec()$models){
 key<-paste(version,cohort,endpoint,kind,sep=':');x<-score_matrix(z,kind,cohort,revision_phenotypes(version));f<-ridge_logistic(x,y,z$patient_id);stopifnot(f$converged);fits[[key]]<-f
 nb<-fed_band_n(nrow(z));eb<-fed_band_n(sum(y));meta[[key]]<-tibble::tibble(site=site,version=version,cohort=cohort,endpoint=endpoint,model=kind,n_lower=nb[1],n_upper=nb[2],events_lower=eb[1],events_upper=eb[2],status='ok')
 for(term in colnames(x)){
 member<-if(all(x[,term]%in%c(0,1)) && term!='(Intercept)')x[,term]==1 else rep(TRUE,nrow(z))
 safe<-any(member)&&sum(y[member])>=aggregate_min_cell(10) && sum(1-y[member])>=aggregate_min_cell(10) && f$covariance[term,term]>0
 rows[[paste(key,term)]]<-tibble::tibble(site=site,version=version,cohort=cohort,endpoint=endpoint,model=kind,term=term,estimate=if(safe)f$beta[term] else NA_real_,std_error=if(safe)sqrt(f$covariance[term,term]) else NA_real_,status=if(safe)'ok' else 'withheld_small_partition_or_absent')
 }
 }
 }
 }
}
fed_write(dplyr::bind_rows(rows),file.path(out,'coefficients.csv'));fed_write(dplyr::bind_rows(meta),file.path(out,'model_metadata.csv'));fed_json(list(site=site,spec_hash=fed_hash(score_spec()),specification=score_spec()),file.path(out,'manifest.json'));saveRDS(fits,file.path(home,'local/score_fits.rds'))
message(site,' severity comparator training complete')
