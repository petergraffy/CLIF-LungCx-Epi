source('diagnostics/revision20261022/lib.R')
site<-commandArgs(trailingOnly=TRUE)[1];home<-file.path(rev_root(),'sites',site);d<-readRDS(file.path(home,'local/review.rds'))$data
# Fraction coefficients are scaled to 24 positive hours (one third of window).
for(v in c('hypox_fraction','hyper_fraction','mixed_fraction'))d[[paste0(v,'_24h')]]<-d[[v]]*3
terms<-paste0(c('hypox_fraction','hyper_fraction','mixed_fraction'),'_24h');rows<-list();wald<-list()
for(endpoint in c('death','hospice','death_or_hospice')){
 d$event<-as.integer(d$time_from_landmark<=28 & if(endpoint=='death_or_hospice')d$event_type%in%c('death','hospice') else d$event_type==endpoint)
 covars<-c(model_spec()$core,paste0('phenotype_P',2:4),terms);covars<-covars[vapply(d[covars],function(x)length(unique(x))>1,logical(1))]
 f<-as.formula(paste('survival::Surv(followup,event)~',paste(c(covars,'strata(info_status)','cluster(patient_id)'),collapse='+')))
 warn<-character();fit<-tryCatch(withCallingHandlers(survival::coxph(f,data=d),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')}),error=identity)
 good<-inherits(fit,'coxph') && !length(warn) && all(is.finite(coef(fit)))
 safe_all<-TRUE
 for(term in terms){
 pos<-d[[term]]>0;safe<-all(vapply(c(FALSE,TRUE),function(v)sum(d$event[pos==v])>=10 && sum(1-d$event[pos==v])>=10,logical(1)));safe_all<-safe_all && safe
 ok<-good && safe && term%in%names(coef(fit));rows[[paste(endpoint,term)]]<-tibble::tibble(site=site,endpoint=endpoint,term=term,estimate=if(ok)unname(coef(fit)[term]) else NA_real_,std_error=if(ok)sqrt(vcov(fit)[term,term]) else NA_real_,status=if(ok)'ok' else if(!safe)'withheld_small_partition' else 'nonestimable',effect='cause-specific log HR per additional 24h evidence, conditional on other channels, phenotype and core covariates')
 }
 good<-good && safe_all && all(terms%in%names(coef(fit)));stat<-if(good)tryCatch(drop(t(coef(fit)[terms])%*%solve(vcov(fit)[terms,terms],coef(fit)[terms])),error=function(e)NA_real_) else NA_real_
 wald[[endpoint]]<-tibble::tibble(site=site,endpoint=endpoint,df=3,robust_wald_statistic=stat,p_value=if(is.finite(stat))pchisq(stat,3,lower.tail=FALSE) else NA_real_,status=if(is.finite(stat))'ok' else 'withheld_or_nonestimable')
}
fed_write(dplyr::bind_rows(rows),file.path(home,'export/arf_associations.csv'));fed_write(dplyr::bind_rows(wald),file.path(home,'export/arf_joint_tests.csv'))
