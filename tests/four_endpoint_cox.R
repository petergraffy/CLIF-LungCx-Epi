source('diagnostics/federated/model_lib.R');source('diagnostics/federated/cox_lib.R')
edge<-data.frame(time_from_landmark=c(1,7,28,28.01,40),event_type=c('death','hospice','live_discharge','live_discharge','death'))
stopifnot(identical(cox_event(edge,'death'),c(1L,0L,0L,0L,0L)),identical(cox_event(edge,'hospice'),c(0L,1L,0L,0L,0L)),identical(cox_event(edge,'death_or_hospice'),c(1L,1L,0L,0L,0L)),identical(cox_event(edge,'live_discharge'),c(0L,0L,1L,0L,0L)))
set.seed(410);n<-600L
d<-tibble::tibble(patient_id=rep(paste0('patient',1:300),each=2),hospitalization_id=paste0('stay',1:n),phenotype=rep(c('P1','P2','P3','P4','I'),length.out=n),age_years=rnorm(n,70,8),sex_category=sample(c('male','female'),n,TRUE),race_category=sample(c('white','black','other'),n,TRUE),admission_year=sample(2018:2025,n,TRUE),charlson_noncancer=sample(0:5,n,TRUE),advanced_cancer=sample(0:1,n,TRUE),code_at_admission=sample(c('full','dnr','unknown'),n,TRUE),sofa_total=sample(c(NA,0:15),n,TRUE),time_from_landmark=rexp(n,1/15)+.1,event_type=sample(c('death','hospice','live_discharge'),n,TRUE,c(.3,.2,.5)))
x<-fed_cox_associations(d,paste0('P',1:4),'synthetic','full','frozen_test')
stopifnot(nrow(x$metadata)==12,setequal(x$metadata$endpoint,cox_endpoints()),all(x$metadata$status=='ok'),all(x$metadata$n_lower[x$metadata$model=='core']==n),all(x$data$followup<=28),all(x$data$info_status[x$data$phenotype=='I']=='insufficient'))
for(endpoint in cox_endpoints()){
 z<-x$data;z$event<-cox_event(z,endpoint)
 formula<-as.formula(paste('survival::Surv(followup,event)~',paste(c(paste0('phenotype_P',2:4),model_spec()$core,'strata(info_status)','cluster(patient_id)'),collapse='+')))
 expected<-survival::coxph(formula,data=z,ties='efron',x=TRUE,y=TRUE,model=FALSE)
 actual<-x$fits[[paste0('core_',endpoint)]]
 stopifnot(max(abs(coef(expected)-coef(actual)))<1e-12,max(abs(vcov(expected)-vcov(actual)))<1e-12,actual$nevent==sum(z$event),actual$n==n)
}
# Diagnostics must reconstruct live-discharge coding, including day28 censoring.
ph<-supplement_model_diagnostics(list(data=x$data,cox=x$fits),x$metadata)
stopifnot(all(ph$table$status=='computed_original_coefficients_reproduced'),nrow(ph$table)==12)
piece<-supplement_piecewise_associations(x$data)
for(endpoint in cox_endpoints()){
 q<-piece$metadata[piece$metadata$endpoint==endpoint,]
 stopifnot(sum(q$events_n)==sum(cox_event(x$data,endpoint)),q$candidate_n[q$interval_start_day==7]==sum(x$data$followup>7))
}
cat('PASS: four cause-specific endpoints, competing exits/day28 boundary, common cohort, clustered Cox agreement, PH reconstruction and piecewise event conservation\n')
