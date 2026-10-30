source("diagnostics/federated/model_lib.R")
args<-commandArgs(trailingOnly=TRUE);site<-args[1];version<-args[2];root<-fed_root();home<-file.path(root,"sites",site)
consensus<-file.path(root,"coordinator",version);manifest<-jsonlite::fromJSON(file.path(consensus,"frozen_manifest.json"));proto<-readRDS(file.path(consensus,"prototypes.rds"))
stopifnot(manifest$spec_hash==fed_hash(fed_spec()))
out<-file.path(home,"export","outcomes",version);private<-file.path(home,"local",version)
if(file.exists(file.path(out,"manifest.json")))stop("Local outcome pass already exists")
dir.create(private,recursive=TRUE,showWarnings=FALSE)
a<-readRDS(file.path(home,"local","discovery.rds"));d<-a$data
features_cache<-file.path(home,"local","patient_features.rds")
if(file.exists(features_cache))x<-readRDS(features_cache) else {x<-fed_patient_features(a$panel,d$hospitalization_id);saveRDS(x,features_cache)}
stopifnot(identical(rownames(x),d$hospitalization_id))
assigned<-fed_assign(x,proto,72-d$hours_UNKNOWN);d<-dplyr::bind_cols(d,assigned)
stopifnot(nrow(d)==nrow(a$data),!anyNA(d$phenotype),all(d$time_from_landmark>0))
d$info_status<-ifelse(d$phenotype=="I","insufficient","classified")
d$followup<-pmin(28,d$time_from_landmark)
base<-base_matrix(d);for(nm in setdiff(colnames(base),"(Intercept)"))d[[nm]]<-base[,nm]
phenotypes<-rownames(proto);pterms<-paste0("phenotype_",phenotypes[-1]);for(j in seq_along(pterms))d[[pterms[j]]]<-as.numeric(d$phenotype==phenotypes[j+1])
composition<-d |> dplyr::count(phenotype,assignment_status,name="n")
bands<-t(vapply(composition$n,fed_band_n,numeric(2)));composition$n_lower<-bands[,1];composition$n_upper<-bands[,2]
fed_write(composition |> dplyr::mutate(site=site,version=version),file.path(out,"assignment_counts.csv"))
assoc<-list();covariances<-list();metadata<-list();fits<-list()
for(spec in c("unadjusted","core","measured_sofa"))for(endpoint in c("death","hospice","death_or_hospice")){
 z<-if(spec=="measured_sofa")d[!is.na(d$sofa_total),] else d
 z$event<-as.integer(z$time_from_landmark<=28 & if(endpoint=="death_or_hospice")z$event_type%in%c("death","hospice") else z$event_type==endpoint)
 covars<-if(spec=="unadjusted")character() else model_spec()$core
 if(spec=="measured_sofa")covars<-c(covars,"sofa_total")
 covars<-covars[vapply(z[covars],function(q)length(unique(q))>1,logical(1))]
 terms<-pterms[vapply(z[pterms],function(q)length(unique(q))>1,logical(1))]
 formula<-as.formula(paste("survival::Surv(followup,event)~",paste(c(terms,covars,"strata(info_status)","cluster(patient_id)"),collapse="+")))
 warnings<-character();fit<-tryCatch(withCallingHandlers(survival::coxph(formula,data=z,ties="efron",x=TRUE,y=TRUE,model=FALSE),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=identity)
 good<-inherits(fit,"coxph") && !length(warnings) && all(is.finite(coef(fit)))
 key<-paste(spec,endpoint,sep="_");fits[[key]]<-fit
 nb<-fed_band_n(nrow(z));eb<-fed_band_n(sum(z$event));metadata[[key]]<-tibble::tibble(site=site,version=version,model=spec,endpoint=endpoint,n_lower=nb[1],n_upper=nb[2],events_lower=eb[1],events_upper=eb[2],status=if(good)"ok" else if(inherits(fit,"error"))conditionMessage(fit) else paste(warnings,collapse=";"),formula=paste(deparse(formula),collapse=" "))
 safe_terms<-character()
 for(j in seq_along(pterms)){
 ph<-phenotypes[j+1];groups<-c("P1",ph);safe<-all(vapply(groups,function(g){zz<-z[z$phenotype==g,];sum(zz$event)>=aggregate_min_cell(10) && sum(!zz$event)>=aggregate_min_cell(10)},logical(1)))
 term<-pterms[j];estimable<-good && term%in%names(coef(fit)) && safe
 status<-if(!good)"model_nonestimable" else if(!safe)"withheld_small_event_partition" else if(!term%in%names(coef(fit)))"phenotype_absent" else "ok"
 assoc[[paste(key,ph)]]<-tibble::tibble(site=site,version=version,consensus_id=manifest$consensus_id,model=spec,endpoint=endpoint,phenotype=ph,reference="P1",estimate=if(estimable)unname(coef(fit)[term]) else NA_real_,std_error=if(estimable)sqrt(vcov(fit)[term,term]) else NA_real_,status=status,effect_type="28day_cause_specific_log_HR")
 if(estimable)safe_terms<-c(safe_terms,term)
 }
 if(length(safe_terms)){v<-as.data.frame(as.table(vcov(fit)[safe_terms,safe_terms,drop=FALSE]));names(v)<-c("term1","term2","covariance");covariances[[key]]<-v |> dplyr::mutate(site=site,version=version,model=spec,endpoint=endpoint)}
}
fed_write(dplyr::bind_rows(assoc),file.path(out,"associations.csv"));fed_write(dplyr::bind_rows(metadata),file.path(out,"model_metadata.csv"))
fed_write(dplyr::bind_rows(covariances),file.path(out,"association_covariances.csv"))
# Local Aalen-Johansen curves. Early/rare event partitions are masked jointly.
source('diagnostics/manuscript_figures/cif_lib.R')
cif<-released_cif_grid(d,days=c(0,1,3,7,14,28))
fed_write(cif |> dplyr::mutate(site=site,version=version),file.path(out,"cumulative_incidence.csv"))
# Predictive models are learned locally; exported coefficients have no IDs.
y<-as.integer(d$time_from_landmark<=28 & d$event_type%in%c("death","hospice"));prediction<-list();predfits<-list()
for(kind in c("baseline","plus_phenotype")){
 xx<-model_matrix(d,phenotypes,plus=kind=="plus_phenotype");f<-ridge_logistic(xx,y,d$patient_id);stopifnot(f$converged);predfits[[kind]]<-f
 for(term in colnames(xx)){
 member<-if(grepl("^phenotype_",term))d$phenotype==sub("^phenotype_","",term) else if(all(xx[,term]%in%c(0,1)) && term!="(Intercept)")xx[,term]==1 else rep(TRUE,nrow(d))
 safe<-any(member)&&sum(y[member])>=aggregate_min_cell(10) && sum(1-y[member])>=aggregate_min_cell(10) && f$covariance[term,term]>0
 prediction[[paste(kind,term)]]<-tibble::tibble(site=site,version=version,model=kind,term=term,estimate=if(safe)f$beta[term] else NA_real_,std_error=if(safe)sqrt(f$covariance[term,term]) else NA_real_,status=if(safe)"ok" else "withheld_small_partition_or_absent")
 }
}
fed_write(dplyr::bind_rows(prediction),file.path(out,"prediction_coefficients.csv"))
fed_json(list(schema=FED_SCHEMA,site=site,version=version,consensus_id=manifest$consensus_id,consensus_hash=fed_hash(manifest),
 specification=model_spec(),model_spec_hash=fed_hash(model_spec()),candidate_hospitalizations=nrow(d),
 training_role=if(site%in%manifest$training_sites)"training_site" else "held_out_site",source=as.list(tools::md5sum(c("diagnostics/federated/model_lib.R","diagnostics/federated/site_outcomes.R")))),file.path(out,"manifest.json"))
saveRDS(list(data=d,assignment=assigned,cox=fits,predictive=predfits),file.path(private,"outcomes.rds"))
fed_json(list(n=nrow(d),all_retained=TRUE,phenotype_counts=as.list(table(d$phenotype)),assignment_status=as.list(table(d$assignment_status))),file.path(private,"outcome_audit.json"))
message(site," ",version," outcome pass complete; all ",nrow(d)," hospitalizations retained")
