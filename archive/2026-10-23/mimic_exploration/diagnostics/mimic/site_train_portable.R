source('diagnostics/mimic/portable_lib.R')
site<-commandArgs(trailingOnly=TRUE)[1];stopifnot(site%in%portable_spec()$training_sites)
d<-readRDS(file.path('output/run_federated_rehearsal_20261022/sites',site,'local/full/outcomes.rds'))$data
rows<-list();meta<-list()
for(cohort in c('all','measured_sofa'))for(endpoint in c('death_or_hospice','death'))for(kind in if(cohort=='all')c('baseline','phenotype','summaries','summaries_phenotype') else portable_spec()$models){
 z<-if(cohort=='all')d else d[!is.na(d$sofa_total),];y<-as.integer(z$time_from_landmark<=28 & if(endpoint=='death')z$event_type=='death' else z$event_type%in%c('death','hospice'));x<-portable_matrix(z,kind);f<-ridge_logistic(x,y,z$patient_id);stopifnot(f$converged)
 nb<-fed_band_n(nrow(z));eb<-fed_band_n(sum(y));key<-paste(cohort,endpoint,kind);meta[[key]]<-tibble::tibble(site=site,cohort=cohort,endpoint=endpoint,model=kind,n_lower=nb[1],n_upper=nb[2],events_lower=eb[1],events_upper=eb[2])
 for(term in colnames(x)){
 member<-if(all(x[,term]%in%c(0,1))&&term!='(Intercept)')x[,term]==1 else rep(TRUE,nrow(z))
 safe<-sum(y[member])>=10 && sum(1-y[member])>=10 && f$covariance[term,term]>0
 rows[[paste(key,term)]]<-tibble::tibble(site=site,cohort=cohort,endpoint=endpoint,model=kind,term=term,estimate=if(safe)f$beta[term] else NA_real_,std_error=if(safe)sqrt(f$covariance[term,term]) else NA_real_,status=if(safe)'ok' else 'withheld_small_partition_or_absent')
 }
}
out<-file.path(mimic_root(),'training',site,'export');fed_write(dplyr::bind_rows(rows),file.path(out,'coefficients.csv'));fed_write(dplyr::bind_rows(meta),file.path(out,'metadata.csv'));fed_json(list(site=site,spec=portable_spec(),spec_hash=fed_hash(portable_spec())),file.path(out,'manifest.json'))
