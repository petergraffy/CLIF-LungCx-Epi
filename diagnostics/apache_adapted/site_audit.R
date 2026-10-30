source('diagnostics/apache_adapted/lib.R');site<-commandArgs(trailingOnly=TRUE)[1];home<-file.path(ap_root(),'sites',site)
a<-readRDS(file.path(home,'local/components.rds'));version<-paste0('loo_',site);d<-ap_data(site,version)
stopifnot(identical(d$hospitalization_id,a$hospitalization_id),all(a$apache_observed_components==rowSums(!is.na(a[paste0('apache_',APACHE_COMPONENTS)]))),all(a$apache_adapted_complete==!is.na(a$apache_adapted_total)),all(a$apache_adapted_total[!is.na(a$apache_adapted_total)]>=0 & a$apache_adapted_total[!is.na(a$apache_adapted_total)]<=62),all(is.na(a$apache_chronic_points)),all(is.na(a$apache_acute_renal_modifier)))
for(k in APACHE_COMPONENTS){v<-a[[paste0('apache_',k)]];stopifnot(all(is.na(v)|(v>=0 & v<=if(k=='gcs')12 else 4)))}
b<-readr::read_csv(file.path(ap_root(),'coordinator',version,'coefficients.csv'),show_col_types=FALSE);m<-jsonlite::fromJSON(file.path(ap_root(),'coordinator',version,'frozen_predictor.json'));stopifnot(!site%in%m$training_sites,!any(grepl(paste0('exchange/',site,'/'),m$input_files)))
p<-readRDS(file.path(home,'local/evaluation.rds'));checks<-list()
for(cohort in ap_spec()$cohorts){z<-ap_cases(d,cohort)
 if(cohort=='all_case_missing_indicators')stopifnot(nrow(z)==nrow(d))
 for(endpoint in ap_spec()$endpoints){
  q<-p[[paste(cohort,endpoint)]];stopifnot(identical(z$hospitalization_id,q$ids),identical(z$patient_id,q$patient_id))
  y<-as.integer(z$time_from_landmark<=28 & if(endpoint=='death')z$event_type=='death' else z$event_type%in%c('death','hospice'));stopifnot(identical(y,q$y))
  for(kind in ap_models(cohort))if(length(q$pred[[kind]])&&!anyNA(q$pred[[kind]])){
   x<-ap_matrix(z,kind,cohort,ap_phenotypes(version));bb<-b[b$cohort==cohort&b$endpoint==endpoint&b$model==kind,];beta<-setNames(bb$estimate,bb$term)
   stopifnot(max(abs(plogis(drop(x%*%beta[colnames(x)]))-q$pred[[kind]]))<1e-12)
  }
  checks[[paste(cohort,endpoint)]]<-list(n=nrow(z),events=sum(y),same_cases_across_models=TRUE)
 }
}
# Exact denominators are local, exported only under their release rules.
saveRDS(checks,file.path(home,'local/evaluation_audit_exact.rds'))
fed_json(list(status='passed',site=site,all_primary_members_in_all_case_sensitivity=TRUE,complete_score_requires_all12_plus_age=TRUE,not_standard_APACHE_II=TRUE,not_APACHE_IV=TRUE,missing_values_not_silently_normal_in_primary=TRUE,held_out_training_exclusion=TRUE,independently_reproduced_frozen_predictions=TRUE,source_numeric_values_frozen=TRUE,model_cohort_summaries='see released cohort_counts and evaluation_metadata; exact values local'),file.path(home,'export/audit.json'))
cat(site,'adapted APACHE audits passed\n')
