source('diagnostics/revision20261022/score_lib.R')
ap_root<-function()normalizePath(Sys.getenv('CLIF_APACHE_OUTPUT',unset='output/run_APACHE_adapted_20261030_v2'),mustWork=FALSE)
ap_primary<-function()normalizePath(Sys.getenv('CLIF_APACHE_PRIMARY',unset='output/run_federated_rehearsal_20261026'))
ap_sites<-function()strsplit(Sys.getenv('CLIF_APACHE_SITES',unset='UCMC,NU,RUSH'),',',fixed=TRUE)[[1]]
ap_spec<-function()list(score=apache_adapted_spec(),sites=ap_sites(),cohorts=c('complete_hb_proxy','all_case_missing_indicators','matched_complete_SOFA'),endpoints=c('death','death_or_hospice'),horizon_days=28,bootstrap=200,lambda=1,clinical_validation=FALSE,upstream_review_commit='e71e6a36130e4e011c652fd93ce9433d22a86cb6',training='coefficientwise inverse variance from other-site local ridge models; unsupported coefficients explicitly flagged',scope='new score comparison sensitivity; original phenotypes/outcomes unchanged')
ap_models<-function(cohort)c('baseline','phenotype','adapted_only','baseline_adapted','summaries','summaries_phenotype','adapted_summaries','adapted_summaries_phenotype',if(cohort=='matched_complete_SOFA')c('sofa_only','baseline_sofa','sofa_summaries','sofa_summaries_phenotype'))
ap_cases<-function(d,cohort){if(cohort=='complete_hb_proxy')d[d$apache_adapted_complete,] else if(cohort=='matched_complete_SOFA')d[d$apache_adapted_complete & !is.na(d$sofa_total),] else d}
ap_matrix<-function(d,kind,cohort,phenotypes){
 partial<-cohort=='all_case_missing_indicators';score<-if(partial)d$apache_partial_sum else d$apache_adapted_total
 if(kind=='adapted_only')x<-cbind(`(Intercept)`=1,adapted_APS_age=score/10) else if(kind=='sofa_only')return(cbind(`(Intercept)`=1,modified_SOFA=d$sofa_total/5)) else {
  x<-base_matrix(d)
  if(grepl('adapted',kind))x<-cbind(x,adapted_APS_age=score/10)
  if(grepl('sofa',kind))x<-cbind(x,modified_SOFA=d$sofa_total/5)
  if(grepl('summaries',kind)){
   x<-cbind(x,imv_fraction=d$hours_IMV/72,hfnc_fraction=d$hours_HFNC/72,cpap_fraction=d$hours_CPAP/72,nippv_fraction=d$hours_NIPPV/72,low_o2_fraction=d$hours_LOW_O2/72,arf_positive=as.numeric(d$measured_arf),unknown_support_fraction=d$hours_UNKNOWN/72,unknown_physiology_fraction=d$physiology_unknown_hours/72)
   for(st in setdiff(SUPPORT_STATES,'ROOM_AIR'))x<-cbind(x,setNames(data.frame(as.numeric(d$support_at_landmark==st)),paste0('hour72_',st))) |> as.matrix()
  }
  if(kind=='phenotype'||grepl('_phenotype$',kind))for(ph in c(setdiff(phenotypes,'P1'),'I'))x<-cbind(x,setNames(data.frame(as.numeric(d$phenotype==ph)),paste0('phenotype_',ph))) |> as.matrix()
 }
 if(partial&&grepl('adapted',kind)){
  mi<-as.matrix(is.na(d[paste0('apache_',APACHE_COMPONENTS)]))*1;colnames(mi)<-paste0('missing_',APACHE_COMPONENTS);x<-cbind(x,mi)
 }
 stopifnot(!anyNA(x),all(is.finite(x)));x
}
ap_data<-function(site,version){
 d<-readRDS(file.path(ap_primary(),'sites',site,'local',version,'outcomes.rds'))$data
 a<-readRDS(file.path(ap_root(),'sites',site,'local/components.rds'))
 d<-dplyr::left_join(d,a,by='hospitalization_id');stopifnot(!anyNA(d$apache_adapted_complete));d
}
ap_pairs<-function(cohort)c(list(phenotype_beyond_adapted_and_summaries=c('adapted_summaries_phenotype','adapted_summaries'),phenotype_beyond_adapted=c('adapted_summaries_phenotype','baseline_adapted'),phenotype_beyond_summaries=c('summaries_phenotype','summaries'),adapted_beyond_baseline=c('baseline_adapted','baseline')),if(cohort=='matched_complete_SOFA')list(adapted_vs_SOFA=c('adapted_only','sofa_only'),baseline_adapted_vs_baseline_SOFA=c('baseline_adapted','baseline_sofa'),adapted_summaries_vs_SOFA_summaries=c('adapted_summaries','sofa_summaries'),phenotype_beyond_SOFA_and_summaries=c('sofa_summaries_phenotype','sofa_summaries')))
ap_export<-function(x,path)fed_write(x,path)
ap_phenotypes<-function(version)rownames(readRDS(file.path(ap_primary(),'coordinator',version,'prototypes.rds')))
