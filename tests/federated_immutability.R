source("diagnostics/federated/comparator_lib.R")
root<-fed_root();check<-function(x,label){if(!isTRUE(x))stop(label);cat("PASS:",label,"\n")}
ledger<-list()
for(version in c("full","loo_UCMC","loo_NU","loo_RUSH")){
 path<-file.path(root,"coordinator",version,"frozen_manifest.json");m<-jsonlite::fromJSON(path,simplifyVector=FALSE)
 for(group in m$source_hashes)for(file in names(group))stopifnot(unname(tools::md5sum(file))==group[[file]])
 proto<-readRDS(file.path(root,"coordinator",version,"prototypes.rds"))
 expected<-substr(fed_hash(list(spec=fed_spec(),training_sites=as.character(unlist(m$training_sites)),source=m$source_hashes,prototypes=proto)),1,20)
 check(expected==m$consensus_id,paste("frozen consensus content hash intact:",version))
 simplified<-jsonlite::fromJSON(path)
 for(site in c("UCMC","NU","RUSH")){
 om<-jsonlite::fromJSON(file.path(root,"exchange","outcomes",version,site,"manifest.json"))
 stopifnot(om$consensus_id==m$consensus_id,om$consensus_hash==fed_hash(simplified))
 }
 for(file in c("frozen_manifest.json","prototypes.rds","outcomes/frozen_predictor.json","outcomes/frozen_comparator_predictor.json")){
 p<-file.path(root,"coordinator",version,file);ledger[[paste(version,file)]]<-list(path=p,md5=unname(tools::md5sum(p)))
 }
}
check(TRUE,"all twelve local outcome passes used the unchanged frozen consensus")
# The discovery snapshot records the exact source before an assignment-only
# matrix-shape fix; assignment is not invoked during local discovery.
archive<-file.path(root,"provenance","discovery_source")
if(dir.exists(archive)){
 original<-c("diagnostics/federated/contract.R","diagnostics/federated/site_discover.R")
 h<-tools::md5sum(file.path(archive,basename(original)));names(h)<-original
 for(site in c("UCMC","NU","RUSH")){m<-jsonlite::fromJSON(file.path(root,"exchange","discovery",site,"manifest.json"));stopifnot(fed_hash(as.list(h))==m$method_hash)}
 check(TRUE,"completed discovery source snapshot matches all site manifests")
}
for(site in c("UCMC","NU","RUSH")){
 a<-readr::read_csv(file.path(root,"exchange","transport",site,"metrics.csv"),show_col_types=FALSE)
 b<-readr::read_csv(file.path(root,"exchange","transport",site,"extended_metrics.csv"),show_col_types=FALSE)
 z<-dplyr::inner_join(a,b,by=c("site","version","model","metric"),suffix=c("_old","_extended"))
 check(max(abs(z$estimate_old-z$estimate_extended))<1e-12,paste("harder comparator leaves original held-out results unchanged:",site))
}
fed_json(list(status="passed",ledger=ledger),file.path(root,"reports","frozen_artifact_ledger.json"))
cat("Frozen-source and train/test artifact checks passed\n")
