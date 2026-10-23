source("diagnostics/federated/model_lib.R")
root<-fed_root();expected<-c("assignment_counts.csv","associations.csv","model_metadata.csv","association_covariances.csv","cumulative_incidence.csv","prediction_coefficients.csv","manifest.json")
for(site in c("UCMC","NU","RUSH"))for(version in c("full","loo_UCMC","loo_NU","loo_RUSH")){
 from<-file.path(root,"sites",site,"export","outcomes",version);to<-file.path(root,"exchange","outcomes",version,site)
 stopifnot(setequal(list.files(from),expected));m<-jsonlite::fromJSON(file.path(from,"manifest.json"));stopifnot(m$site==site,m$version==version,m$model_spec_hash==fed_hash(model_spec()))
 dir.create(to,recursive=TRUE,showWarnings=FALSE)
 for(name in expected){p<-file.path(from,name)
 if(grepl("csv$",name)){z<-readr::read_csv(p,show_col_types=FALSE);stopifnot(!any(names(z)%in%c("patient_id","hospitalization_id","t0","event_time","admission_dttm","discharge_dttm")))
 if(name=="associations.csv")stopifnot(all(is.na(z$estimate[z$status!="ok"])))
 if(name=="assignment_counts.csv")stopifnot(all(z$n_upper-z$n_lower==9))
 }
 stopifnot(file.copy(p,file.path(to,name),overwrite=FALSE))
 }
}
fed_json(list(status="passed",no_patient_identifiers_or_timestamps=TRUE,model_partitions_checked=TRUE,
 source_contract=as.list(tools::md5sum("diagnostics/federated/contract.R"))),file.path(root,"exchange","outcome_gate.json"))
