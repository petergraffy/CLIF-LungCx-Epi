source('diagnostics/revision20261022/lib.R')
score_spec<-function()list(horizon=28,cohorts=c('measured_sofa','missing_domain_zero_sensitivity'),endpoints=c('death_or_hospice','death'),models=c('baseline','phenotype','sofa_only','sofa','summaries','summaries_phenotype','sofa_summaries','sofa_summaries_phenotype'),sofa='first ICU 24h modified measured-component SOFA; creatinine only renal; incomplete scores not standard SOFA',missing='complete-case comparison plus explicitly labeled zero-domain/missing-indicator sensitivity; all primary cases retained outside this sensitivity',apache='II/IV not calculable from current extracts; no fabricated comparison',validation='two-site frozen predictor tested on third; patient bootstrap conditional on trained model',lambda=1)
score_matrix<-function(d,kind,cohort,phenotypes=c('P1','P2','P3','P4')){
 domains<-paste0('sofa_',c('resp','cv','coag','liver','renal','cns'))
 score<-if(cohort=='measured_sofa')d$sofa_total else rowSums(d[domains],na.rm=TRUE)
 if(kind=='sofa_only')return(cbind(`(Intercept)`=1,sofa_score=score/5))
 x<-base_matrix(d)
 if(grepl('sofa',kind)){
 x<-cbind(x,sofa_score=score/5)
 if(cohort!='measured_sofa'){mi<-as.matrix(is.na(d[domains]))*1;colnames(mi)<-paste0(domains,'_missing');x<-cbind(x,mi)}
 }
 if(grepl('summaries',kind)){
 x<-cbind(x,imv_fraction=d$hours_IMV/72,hfnc_fraction=d$hours_HFNC/72,cpap_fraction=d$hours_CPAP/72,nippv_fraction=d$hours_NIPPV/72,low_o2_fraction=d$hours_LOW_O2/72,arf_positive=as.numeric(d$measured_arf),unknown_support_fraction=d$hours_UNKNOWN/72,unknown_physiology_fraction=d$physiology_unknown_hours/72)
 for(st in setdiff(SUPPORT_STATES,'ROOM_AIR'))x<-cbind(x,setNames(data.frame(as.numeric(d$support_at_landmark==st)),paste0('hour72_',st))) |> as.matrix()
 }
 if(kind=='phenotype' || grepl('_phenotype$',kind)){
 for(ph in c(setdiff(phenotypes,'P1'),'I'))x<-cbind(x,setNames(data.frame(as.numeric(d$phenotype==ph)),paste0('phenotype_',ph))) |> as.matrix()
 }
 stopifnot(!anyNA(x));x
}
