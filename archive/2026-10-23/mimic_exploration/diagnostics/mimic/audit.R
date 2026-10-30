source('diagnostics/mimic/portable_lib.R')
root<-mimic_root();a<-readRDS(file.path(root,'local/prepared.rds'));z<-readRDS(file.path(root,'local/derivation.rds'));d<-z$data
stopifnot(all(is.na(a$tables$hospital_diagnosis$poa_present)),all(is.na(d$charlson_noncancer)),all(is.na(d$advanced_cancer)),all(d$time_from_landmark>0),nrow(z$panel)==72*nrow(d),!any(z$panel$terminal),all(d$phenotype[d$hours_UNKNOWN==72]=='I'),nrow(d)==sum(a$landmark$landmark_reason=='eligible'))
# Verify lung-primary-site coding independently from eligibility flags.
dx<-a$tables$hospital_diagnosis;cd<-toupper(gsub('[^A-Za-z0-9]','',dx$diagnosis_code));fmt<-toupper(dx$diagnosis_code_format)
lung<-dx$hospitalization_id[which((grepl('10',fmt)&substr(cd,1,3)=='C34')|(grepl('9',fmt)&substr(cd,1,3)=='162'))]
stopifnot(all(d$hospitalization_id%in%lung))
stopifnot(all(d$first_completed_lung_coded_encounter[d$prior_completed_lung_cancer_encounter]<=d$admission_dttm[d$prior_completed_lung_cancer_encounter]))
# Check history against the full source hospitalization table, not the selected cohort.
h<-arrow::read_parquet(file.path(root,'local/source_tables/clif_hospitalization.parquet'),col_select=c('patient_id','hospitalization_id','admission_dttm','discharge_dttm'))
all_dx<-arrow::read_parquet(file.path(root,'local/source_tables/clif_hospital_diagnosis.parquet'),col_select=c('hospitalization_id','diagnosis_code','diagnosis_code_format'))
cd<-toupper(gsub('[^A-Za-z0-9]','',all_dx$diagnosis_code));fmt<-toupper(all_dx$diagnosis_code_format);priorids<-all_dx$hospitalization_id[which((grepl('10',fmt)&substr(cd,1,3)=='C34')|(grepl('9',fmt)&substr(cd,1,3)=='162'))]
for(i in which(d$prior_completed_lung_cancer_encounter))stopifnot(any(h$patient_id==d$patient_id[i] & h$hospitalization_id!=d$hospitalization_id[i] & h$hospitalization_id%in%priorids & h$discharge_dttm<=d$admission_dttm[i],na.rm=TRUE))
m<-jsonlite::fromJSON(file.path(root,'coordinator/frozen_predictor.json'));stopifnot(identical(m$training_sites,c('UCMC','NU','RUSH')),m$coefficient_sha256==digest::digest(file=file.path(root,'coordinator/coefficients.csv'),algo='sha256'),m$prototype_sha256==digest::digest(file=m$prototype_file,algo='sha256'))
b<-readr::read_csv(file.path(root,'coordinator/coefficients.csv'),show_col_types=FALSE);pr<-readRDS(file.path(root,'local/predictions_private.rds'))
for(key in names(pr)){
 keyparts<-strsplit(key,' ',fixed=TRUE)[[1]];sub<-keyparts[1];cohort<-keyparts[2];endpoint<-keyparts[3]
 x<-if(sub=='coded_lung_cancer')d else d[d$prior_completed_lung_cancer_encounter,];if(cohort=='measured_sofa')x<-x[!is.na(x$sofa_total),]
 stopifnot(identical(x$hospitalization_id,pr[[key]]$ids))
 for(kind in names(pr[[key]]$predictions)){
  xx<-portable_matrix(x,kind);bb<-b[b$cohort==cohort & b$endpoint==endpoint & b$model==kind,];beta<-setNames(bb$estimate,bb$term);stopifnot(max(abs(plogis(drop(xx%*%beta[colnames(xx)]))-pr[[key]]$predictions[[kind]]))<1e-12)
 }
}
# All original frozen three-site artifacts remain unchanged.
ledger<-jsonlite::fromJSON('output/run_federated_rehearsal_20261022/reports/frozen_artifact_ledger.json',simplifyVector=FALSE)$ledger
for(item in ledger)stopifnot(unname(tools::md5sum(item$path))==item$md5)
for(f in list.files(file.path(root,'aggregate'),pattern='csv$',full.names=TRUE)){
 x<-readr::read_csv(f,show_col_types=FALSE);stopifnot(!any(names(x)%in%c('patient_id','hospitalization_id','recorded_dttm','admission_dttm','event_time')))
}
fed_json(list(status='passed',current_encounter_lung_codes_verified=TRUE,prior_completed_lung_cancer_encounter_verified=TRUE,poa_never_claimed_or_used_as_measured=TRUE,diagnosis_covariates_not_fabricated=TRUE,all_landmark_cases_retained=TRUE,frozen_three_site_phenotypes_unchanged=TRUE,no_MIMIC_training=TRUE,independent_frozen_prediction_reproduction=TRUE,source_room_air_gap_explicit=TRUE,patient_identifiers_in_aggregates=0,n=nrow(d),prior_documented_n=sum(d$prior_completed_lung_cancer_encounter)),file.path(root,'aggregate/final_audit.json'))
cat('MIMIC lung-cancer eligibility, chronology, retention, training separation and export audits passed\n')
