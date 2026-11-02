source("diagnostics/federated/comparator_lib.R")
site<-commandArgs(trailingOnly=TRUE)[1];root<-fed_root();home<-file.path(root,"sites",site)
for(version in fed_versions()){
 manifest<-jsonlite::fromJSON(file.path(root,"coordinator",version,"frozen_manifest.json"));phenotypes<-paste0("P",1:manifest$k)
 d<-readRDS(file.path(home,"local",version,"outcomes.rds"))$data;y<-as.integer(d$time_from_landmark<=28 & d$event_type%in%c("death","hospice"))
 rows<-list()
 for(kind in c("simple_summaries","summaries_plus_phenotype")){
 x<-comparator_matrix(d,phenotypes,plus=kind=="summaries_plus_phenotype");f<-ridge_logistic(x,y,d$patient_id);stopifnot(f$converged)
 for(term in colnames(x)){
 member<-if(grepl("^phenotype_",term))d$phenotype==sub("^phenotype_","",term) else if(all(x[,term]%in%c(0,1)) && term!="(Intercept)")x[,term]==1 else rep(TRUE,nrow(d))
 safe<-any(member)&&sum(y[member])>=aggregate_min_cell(10) && sum(1-y[member])>=aggregate_min_cell(10) && f$covariance[term,term]>0
 rows[[paste(kind,term)]]<-tibble::tibble(site=site,version=version,model=kind,term=term,estimate=if(safe)f$beta[term] else NA_real_,std_error=if(safe)sqrt(f$covariance[term,term]) else NA_real_,status=if(safe)"ok" else "withheld_small_partition_or_absent")
 }
 }
 out<-file.path(home,"export","comparators",version);fed_write(dplyr::bind_rows(rows),file.path(out,"coefficients.csv"))
 fed_json(list(schema=FED_SCHEMA,site=site,version=version,consensus_id=manifest$consensus_id,comparator_spec=comparator_spec(),spec_hash=fed_hash(comparator_spec())),file.path(out,"manifest.json"))
}
message(site," comparator models complete; frozen clustering unchanged")
