source('diagnostics/mimic/portable_lib.R')
root<-mimic_root();a<-readRDS(file.path(root,'local/derivation.rds'));d<-a$data
base<-portable_matrix(d,'baseline');for(term in setdiff(colnames(base),'(Intercept)'))d[[term]]<-base[,term]
for(ph in paste0('P',2:4))d[[paste0('phenotype_',ph)]]<-as.numeric(d$phenotype==ph)
rows<-list();meta<-list()
for(sub in portable_spec()$subgroups)for(endpoint in c('death','hospice','death_or_hospice')){
 z<-if(sub=='coded_lung_cancer')d else d[d$prior_completed_lung_cancer_encounter,]
 z$event<-as.integer(z$time_from_landmark<=28 & if(endpoint=='death_or_hospice')z$event_type%in%c('death','hospice') else z$event_type==endpoint)
 terms<-c(paste0('phenotype_P',2:4),portable_spec()$baseline);terms<-terms[vapply(z[terms],function(x)length(unique(x))>1,logical(1))]
 formula<-as.formula(paste('survival::Surv(followup,event)~',paste(c(terms,'strata(info_status)','cluster(patient_id)'),collapse='+')))
 warnings<-character();f<-tryCatch(withCallingHandlers(survival::coxph(formula,data=z),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=identity);good<-inherits(f,'coxph')&&!length(warnings)&&all(is.finite(coef(f)))
 for(ph in paste0('P',2:4)){
  safe<-all(vapply(c('P1',ph),function(g){q<-z[z$phenotype==g,];sum(q$event)>=10 && sum(1-q$event)>=10},logical(1)));term<-paste0('phenotype_',ph);ok<-good&&safe&&term%in%names(coef(f))
  rows[[paste(sub,endpoint,ph)]]<-tibble::tibble(site='MIMIC',subgroup=sub,endpoint=endpoint,phenotype=ph,reference='P1',estimate=if(ok)unname(coef(f)[term]) else NA_real_,std_error=if(ok)sqrt(vcov(f)[term,term]) else NA_real_,status=if(!any(z$phenotype=='P1'))'reference_phenotype_absent' else if(ok)'ok' else if(!safe)'withheld_small_event_partition' else 'nonestimable',effect='28-day cause-specific log HR; portable demographics/code-status adjustment; not pooled with POA cohort effects')
 }
 nb<-fed_band_n(nrow(z));eb<-fed_band_n(sum(z$event));meta[[paste(sub,endpoint)]]<-tibble::tibble(subgroup=sub,endpoint=endpoint,n_lower=nb[1],n_upper=nb[2],events_lower=eb[1],events_upper=eb[2],status=if(!any(z$phenotype=='P1'))'reference_phenotype_absent' else if(good)'ok' else 'nonestimable',formula=paste(deparse(formula),collapse=' '))
}
fed_write(dplyr::bind_rows(rows),file.path(root,'aggregate/associations.csv'));fed_write(dplyr::bind_rows(meta),file.path(root,'aggregate/association_metadata.csv'))
