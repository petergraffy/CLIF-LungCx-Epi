source('diagnostics/federated/model_lib.R')
mimic_root<-function()Sys.getenv('CLIF_MIMIC_OUTPUT',unset='output/run_MIMIC_external_20261023')
portable_spec<-function()list(schema='mimic-external-portable-v1',training_sites=c('UCMC','NU','RUSH'),external_site='MIMIC',horizon=28L,cohort='current-encounter C34/162 lung cancer, adult ICU +72h landmark; prior completed cancer encounter subgroup',baseline=c('age10','male','black','other_or_unknown_race','code_nonfull','code_unknown'),excluded='calendar year (shifted); POA Charlson/advanced-cancer proxy (not observed); discharge diagnoses are not baseline predictors',primary_consensus_id=jsonlite::fromJSON('output/run_federated_rehearsal_20261022/coordinator/full/frozen_manifest.json')$consensus_id,models=c('baseline','phenotype','summaries','summaries_phenotype','sofa_only','sofa','sofa_summaries','sofa_summaries_phenotype'),lambda=1L,bootstrap=200L,subgroups=c('coded_lung_cancer','prior_completed_lung_cancer_encounter'),missing='unknown support retained; complete-component modified-SOFA cohort separately labeled; no all-case score imputation in this external comparison',freeze='three-site templates unchanged; three-site portable coefficients frozen before MIMIC prediction evaluation; no MIMIC refit or recalibration')
portable_matrix<-function(d,kind){
 x<-cbind(`(Intercept)`=1,age10=(d$age_years-70)/10,male=as.integer(tolower(d$sex_category)=='male'),black=as.integer(grepl('black|african',tolower(d$race_category))),other_or_unknown_race=as.integer(!grepl('black|african|white',tolower(d$race_category))),code_nonfull=as.integer(!d$code_at_admission%in%c('full','unknown')),code_unknown=as.integer(d$code_at_admission=='unknown'))
 if(kind=='sofa_only')return(cbind(`(Intercept)`=1,sofa_score=d$sofa_total/5))
 if(grepl('sofa',kind))x<-cbind(x,sofa_score=d$sofa_total/5)
 if(grepl('summaries',kind)){
 x<-cbind(x,imv_fraction=d$hours_IMV/72,hfnc_fraction=d$hours_HFNC/72,cpap_fraction=d$hours_CPAP/72,nippv_fraction=d$hours_NIPPV/72,low_o2_fraction=d$hours_LOW_O2/72,arf_positive=as.numeric(d$measured_arf),unknown_support_fraction=d$hours_UNKNOWN/72,unknown_physiology_fraction=d$physiology_unknown_hours/72)
 for(st in setdiff(SUPPORT_STATES,'ROOM_AIR'))x<-cbind(x,setNames(data.frame(as.numeric(d$support_at_landmark==st)),paste0('hour72_',st))) |> as.matrix()
 }
 if(kind=='phenotype'||grepl('_phenotype$',kind))for(ph in c('P2','P3','P4','I'))x<-cbind(x,setNames(data.frame(as.numeric(d$phenotype==ph)),paste0('phenotype_',ph))) |> as.matrix()
 stopifnot(!anyNA(x));x
}
