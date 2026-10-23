source('diagnostics/revision20261022/score_lib.R')
stopifnot(auc_rank(c(0,0,1,1),c(.1,.2,.8,.9))==1)
d<-data.frame(age_years=c(60,70,80),sex_category=c('Male','Female','Male'),race_category=c('White','Black','Other'),admission_year=c(2019,2020,2021),charlson_noncancer=c(1,2,3),advanced_cancer=c(0,1,0),code_at_admission=c('full','dni','unknown'),sofa_total=c(6,NA,0),sofa_resp=c(2,NA,0),sofa_cv=c(1,0,0),sofa_coag=c(0,1,0),sofa_liver=c(0,1,0),sofa_renal=c(1,1,0),sofa_cns=c(2,1,0),hours_IMV=c(0,72,0),hours_HFNC=0,hours_CPAP=0,hours_NIPPV=0,hours_LOW_O2=c(72,0,0),hours_UNKNOWN=c(0,0,72),measured_arf=c(TRUE,TRUE,FALSE),physiology_unknown_hours=c(5,10,72),support_at_landmark=c('LOW_O2','IMV','UNKNOWN'),phenotype=c('P1','P4','I'))
x<-score_matrix(d,'sofa_summaries_phenotype','missing_domain_zero_sensitivity')
stopifnot(nrow(x)==3,!anyNA(x),x[2,'sofa_resp_missing']==1,x[3,'hour72_UNKNOWN']==1,x[3,'phenotype_I']==1,x[2,'imv_fraction']==1)
# Neither outcomes nor patient IDs are allowed into any predictive design.
d$event_type<-c('death','hospice','live_discharge');d$time_from_landmark<-c(1,2,3)
stopifnot(identical(x,score_matrix(d,'sofa_summaries_phenotype','missing_domain_zero_sensitivity')))
for(kind in score_spec()$models){xx<-score_matrix(d[!is.na(d$sofa_total),],kind,'measured_sofa');stopifnot(nrow(xx)==2,!anyNA(xx),!any(grepl('event|patient|time_from',colnames(xx))))}
stopifnot(all(is.na(d$sofa_total)==c(FALSE,TRUE,FALSE)))
# No incomplete score is labeled APACHE, and no validated APACHE performance is invented.
stopifnot(grepl('not calculable',score_spec()$apache))
cat('Revision design contracts passed\n')
