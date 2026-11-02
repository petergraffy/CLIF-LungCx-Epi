# Shared association-only fitter; predictors and phenotype learning are separate.
# Callers source model_lib.R first.
fed_cox_associations<-function(d,phenotypes,site,version,consensus_id,keep_x=TRUE){
 stopifnot(identical(phenotypes[1],'P1'),!anyNA(d$phenotype),all(d$time_from_landmark>0))
 d$info_status<-ifelse(d$phenotype=='I','insufficient','classified')
 d$followup<-pmin(28,d$time_from_landmark)
 base<-base_matrix(d);for(nm in setdiff(colnames(base),'(Intercept)'))d[[nm]]<-base[,nm]
 pterms<-paste0('phenotype_',phenotypes[-1]);for(j in seq_along(pterms))d[[pterms[j]]]<-as.numeric(d$phenotype==phenotypes[j+1])
assoc<-list();covariances<-list();metadata<-list();fits<-list()
for(spec in c("unadjusted","core","measured_sofa"))for(endpoint in cox_endpoints()){
 z<-if(spec=="measured_sofa")d[!is.na(d$sofa_total),] else d
 z$event<-cox_event(z,endpoint)
 covars<-if(spec=="unadjusted")character() else model_spec()$core
 if(spec=="measured_sofa")covars<-c(covars,"sofa_total")
 covars<-covars[vapply(z[covars],function(q)length(unique(q))>1,logical(1))]
 terms<-pterms[vapply(z[pterms],function(q)length(unique(q))>1,logical(1))]
 formula<-as.formula(paste("survival::Surv(followup,event)~",paste(c(terms,covars,"strata(info_status)","cluster(patient_id)"),collapse="+")))
 environment(formula)<-.GlobalEnv
 warnings<-character();fit<-tryCatch(withCallingHandlers(survival::coxph(formula,data=z,ties="efron",x=keep_x,y=TRUE,model=FALSE),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=identity)
 good<-inherits(fit,"coxph") && !length(warnings) && all(is.finite(coef(fit)))
 key<-paste(spec,endpoint,sep="_");fits[[key]]<-fit
 nb<-fed_band_n(nrow(z));eb<-fed_band_n(sum(z$event));metadata[[key]]<-tibble::tibble(site=site,version=version,model=spec,endpoint=endpoint,n_lower=nb[1],n_upper=nb[2],events_lower=eb[1],events_upper=eb[2],status=if(good)"ok" else if(inherits(fit,"error"))conditionMessage(fit) else paste(warnings,collapse=";"),formula=paste(deparse(formula),collapse=" "))
 safe_terms<-character()
 for(j in seq_along(pterms)){
 ph<-phenotypes[j+1];groups<-c("P1",ph);safe<-all(vapply(groups,function(g){zz<-z[z$phenotype==g,];sum(zz$event)>=aggregate_min_cell(10) && sum(!zz$event)>=aggregate_min_cell(10)},logical(1)))
 term<-pterms[j];estimable<-good && term%in%names(coef(fit)) && safe
 status<-if(!good)"model_nonestimable" else if(!safe)"withheld_small_event_partition" else if(!term%in%names(coef(fit)))"phenotype_absent" else "ok"
 assoc[[paste(key,ph)]]<-tibble::tibble(site=site,version=version,consensus_id=consensus_id,model=spec,endpoint=endpoint,phenotype=ph,reference="P1",estimate=if(estimable)unname(coef(fit)[term]) else NA_real_,std_error=if(estimable)sqrt(vcov(fit)[term,term]) else NA_real_,status=status,effect_type="28day_cause_specific_log_HR")
 if(estimable)safe_terms<-c(safe_terms,term)
 }
 if(length(safe_terms)){v<-as.data.frame(as.table(vcov(fit)[safe_terms,safe_terms,drop=FALSE]));names(v)<-c("term1","term2","covariance");covariances[[key]]<-v |> dplyr::mutate(site=site,version=version,model=spec,endpoint=endpoint)}
}

 list(data=d,associations=dplyr::bind_rows(assoc),metadata=dplyr::bind_rows(metadata),covariances=dplyr::bind_rows(covariances),fits=fits)
}
