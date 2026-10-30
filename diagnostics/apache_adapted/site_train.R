source('diagnostics/apache_adapted/lib.R');site<-commandArgs(trailingOnly=TRUE)[1]
home<-file.path(ap_root(),'sites',site);rows<-list();metadata<-list();exact<-list();fits<-list()
for(version in paste0('loo_',ap_sites())){
 d<-ap_data(site,version)
 for(cohort in ap_spec()$cohorts){z<-ap_cases(d,cohort)
 for(endpoint in ap_spec()$endpoints){
 y<-as.integer(z$time_from_landmark<=28 & if(endpoint=='death')z$event_type=='death' else z$event_type%in%c('death','hospice'))
 for(kind in ap_models(cohort)){
 key<-paste(version,cohort,endpoint,kind,sep=':')
 if(!nrow(z)){template<-d[1,];template$apache_adapted_total<-0;template$apache_partial_sum<-0;template$sofa_total<-0;x<-ap_matrix(template,kind,cohort,ap_phenotypes(version))[0,,drop=FALSE]} else x<-ap_matrix(z,kind,cohort,ap_phenotypes(version))
 f<-if(length(y)>=4&&length(unique(y))==2)tryCatch(ridge_logistic(x,y,z$patient_id),error=identity) else NULL
 good<-is.list(f)&&!inherits(f,'error')&&isTRUE(f$converged);fits[[key]]<-f
 nb<-fed_band_n(nrow(z));eb<-fed_band_n(sum(y));pb<-fed_band_n(dplyr::n_distinct(z$patient_id))
 metadata[[key]]<-tibble::tibble(site=site,version=version,cohort=cohort,endpoint=endpoint,model=kind,n_lower=nb[1],n_upper=nb[2],patients_lower=pb[1],patients_upper=pb[2],events_lower=eb[1],events_upper=eb[2],status=if(good)'ok' else'nonestimable_or_insufficient_training_endpoint',design=paste(colnames(x),collapse=' + '))
 exact[[key]]<-tibble::tibble(version=version,cohort=cohort,endpoint=endpoint,model=kind,n=nrow(z),patients=dplyr::n_distinct(z$patient_id),events=sum(y),status=if(good)'ok' else'nonestimable')
 for(term in colnames(x)){
 member<-if(nrow(x)&&all(x[,term]%in%c(0,1))&&term!='(Intercept)')x[,term]==1 else rep(TRUE,nrow(z))
 safe<-good&&any(member)&&sum(y[member])>=aggregate_min_cell(10)&&sum(1-y[member])>=aggregate_min_cell(10)&&f$covariance[term,term]>0
 rows[[paste(key,term)]]<-tibble::tibble(site=site,version=version,cohort=cohort,endpoint=endpoint,model=kind,term=term,estimate=if(safe)f$beta[term] else NA_real_,std_error=if(safe)sqrt(f$covariance[term,term]) else NA_real_,status=if(safe)'ok' else if(!good)'fit_unavailable' else'withheld_small_partition_or_absent')
 }
 }
 }
 }
}
ap_export(dplyr::bind_rows(rows),file.path(home,'export/coefficients.csv'));ap_export(dplyr::bind_rows(metadata),file.path(home,'export/model_metadata.csv'))
readr::write_csv(dplyr::bind_rows(exact),file.path(home,'local/model_denominators_exact.csv'));saveRDS(fits,file.path(home,'local/training_fits.rds'))
fed_json(list(site=site,specification=ap_spec(),spec_hash=fed_hash(ap_spec()),phenotype_learning='frozen original held-out templates; not relearned by this score extension'),file.path(home,'export/training_manifest.json'))
cat(site,'adapted-score training complete\n')
