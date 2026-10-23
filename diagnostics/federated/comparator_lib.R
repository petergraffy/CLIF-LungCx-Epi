source("diagnostics/federated/model_lib.R")
comparator_spec<-function()list(models=c("simple_summaries","summaries_plus_phenotype"),
 summaries=c("imv_fraction","hfnc_fraction","niv_fraction","low_o2_fraction","arf_positive","unknown_support_fraction","unknown_physiology_fraction"),
 purpose="exploratory harder comparator added during rehearsal; clustering and assignment remain frozen",ridge_lambda=1)
comparator_matrix<-function(d,phenotypes,plus=FALSE){
 x<-model_matrix(d,phenotypes,plus=plus)
 cbind(x,imv_fraction=d$hours_IMV/72,hfnc_fraction=d$hours_HFNC/72,niv_fraction=(d$hours_CPAP+d$hours_NIPPV)/72,
 low_o2_fraction=d$hours_LOW_O2/72,arf_positive=as.numeric(d$measured_arf),unknown_support_fraction=d$hours_UNKNOWN/72,
 unknown_physiology_fraction=d$physiology_unknown_hours/72)
}
