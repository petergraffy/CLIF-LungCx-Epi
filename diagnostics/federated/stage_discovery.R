source("diagnostics/federated/contract.R")
root<-fed_root()
for(site in c("UCMC","NU","RUSH")){
 from<-file.path(root,"sites",site,"export","discovery");to<-file.path(root,"exchange","discovery",site)
 expected<-c("fingerprints.csv","clusters.csv","silhouette.csv","manifest.json")
 stopifnot(setequal(list.files(from),expected));dir.create(to,recursive=TRUE,showWarnings=FALSE)
 for(name in expected){p<-file.path(from,name);if(grepl("csv$",name)){
 x<-readr::read_csv(p,show_col_types=FALSE);stopifnot(!any(names(x)%in%c("patient_id","hospitalization_id","event_type","death_in_hosp","hospice_discharge","event_time")))
 };stopifnot(file.copy(p,file.path(to,name),overwrite=FALSE))}
}
fed_json(list(status="passed",only_allowlisted_aggregate_files=TRUE,patient_rows=0,site_dirs=c("UCMC","NU","RUSH")),file.path(root,"exchange","discovery_gate.json"))
