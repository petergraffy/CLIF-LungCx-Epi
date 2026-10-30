source("diagnostics/federated/comparator_lib.R")
site<-commandArgs(trailingOnly=TRUE)[1];root<-fed_root();home<-file.path(root,"sites",site)
a<-readRDS(file.path(home,"local","discovery.rds"));checks<-list()
for(version in c("full","loo_UCMC","loo_NU","loo_RUSH")){
 d<-readRDS(file.path(home,"local",version,"outcomes.rds"))$data
 stopifnot(setequal(d$hospitalization_id,a$data$hospitalization_id),nrow(d)==nrow(a$data),!anyNA(d$phenotype),all(d$time_from_landmark>0))
 insufficient<-72-d$hours_UNKNOWN<12;stopifnot(all(d$phenotype[insufficient]=="I"))
 stopifnot(all(d$phenotype[!insufficient]!="I"),all(d$assignment_status[!insufficient]=="classified"),!any(d$assignment_status%in%c('ambiguous','poor_match')),
           all(d$phenotype[!insufficient]==d$closest_phenotype[!insufficient]))
 m<-jsonlite::fromJSON(file.path(root,"coordinator",version,"frozen_manifest.json"));stopifnot(all(d$phenotype%in%c(paste0("P",1:m$k),"I")))
 checks[[version]]<-list(cohort_retained=TRUE,all_assigned=TRUE,insufficient_information_retained=TRUE,no_prelandmark_outcomes=TRUE)
}
# Independent recomputation ensures the transported predictions use only the
# frozen two-site coefficients, with no held-out fitting or recalibration.
version<-paste0("loo_",site);saved<-readRDS(file.path(home,"local",version,"transport_extended_private.rds"));d<-readRDS(file.path(home,"local",version,"outcomes.rds"))$data
m<-jsonlite::fromJSON(file.path(root,"coordinator",version,"outcomes","frozen_predictor.json"));c<-jsonlite::fromJSON(file.path(root,"coordinator",version,"outcomes","frozen_comparator_predictor.json"))
stopifnot(!site%in%m$training_sites,!site%in%c$training_sites)
coef<-c(m$coefficients,c$coefficients);ph<-paste0("P",1:jsonlite::fromJSON(file.path(root,"coordinator",version,"frozen_manifest.json"))$k)
for(kind in names(coef)){
 x<-if(kind%in%c("baseline","plus_phenotype"))model_matrix(d,ph,kind=="plus_phenotype") else comparator_matrix(d,ph,kind=="summaries_plus_phenotype")
 b<-unlist(coef[[kind]]);stopifnot(max(abs(saved$predictions[[kind]]-plogis(drop(x%*%b[colnames(x)]))))<1e-12)
}
fed_json(list(status="passed",site=site,candidate_hospitalizations=nrow(a$data),measured_sofa_hospitalizations=sum(!is.na(a$data$sofa_total)),support_unknown_hours=sum(a$panel$support=="UNKNOWN"),total_panel_hours=nrow(a$panel),versions=checks,
 held_out_predictions_independently_recomputed=TRUE,held_out_outcomes_not_used_for_fitting=TRUE,
 source_contract=as.list(tools::md5sum("diagnostics/federated/contract.R"))),file.path(home,"export","quality","audit.json"))
message(site," local cohort and prediction audits passed")
