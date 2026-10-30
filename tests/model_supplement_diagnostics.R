source('code/lib/load.R')
set.seed(72);n<-160L
x<-tibble::tibble(patient_id=paste0('p',seq_len(n)),phenotype=rep(paste0('P',1:4),each=40),time_from_landmark=rep(c(2,7,8,14,28,40,4,20),20),event_type=rep(c('death','hospice','live_discharge','death','hospice','live_discharge','live_discharge','death'),20),info_status='classified',age10=rnorm(n))
x$followup<-pmin(x$time_from_landmark,28)
for(p in paste0('P',2:4))x[[paste0('phenotype_',p)]]<-as.numeric(x$phenotype==p)
for(v in c('male','black','other_or_unknown_race','year','charlson5','advanced','code_nonfull','code_unknown'))x[[v]]<-0
x$event<-as.integer(x$time_from_landmark<=28 & x$event_type%in%c('death','hospice'))
f<-survival::coxph(survival::Surv(followup,event)~phenotype_P2+phenotype_P3+phenotype_P4+age10+strata(info_status)+cluster(patient_id),data=x,x=FALSE,model=FALSE)
a<-supplement_model_diagnostics(list(data=x,cox=list(core_death_or_hospice=f)),tibble::tibble(model='core',endpoint='death_or_hospice',status='ok'))
stopifnot(a$table$status=='computed_original_coefficients_reproduced',a$table$coefficient_max_abs_difference==0,nrow(a$term_tests)>0,all(is.finite(a$term_tests$p)))
p<-supplement_piecewise_associations(x,3)
m<-p$metadata[p$metadata$endpoint=='death_or_hospice',]
early<-sum(x$event & x$time_from_landmark<=7);late<-sum(x$event & x$time_from_landmark>7)
stopifnot(m$events_n[m$interval_start_day==0]==early,m$events_n[m$interval_start_day==7]==late,sum(m$events_n)==sum(x$event),m$candidate_n[m$interval_start_day==7]==sum(x$followup>7),all(p$associations$analysis=='exploratory_posthoc_0to7_7to28_days_after_landmark'))
cat('PASS: PH reconstruction and day7 boundary/risk-set/event conservation\n')
