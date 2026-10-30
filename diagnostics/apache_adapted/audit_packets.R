source('diagnostics/apache_adapted/lib.R')
root<-ap_root();counts<-readr::read_csv(file.path(root,'reports/tables/cohort_counts.csv'),show_col_types=FALSE)
files_n<-list()
for(site in ap_sites()){
 p<-file.path(root,'site_packets',site,'share');m<-jsonlite::fromJSON(file.path(p,'site_return_manifest.json'))
 stopifnot(m$apache_adapted_extension_complete,m$apache_not_standard_II_or_IV)
 ledger<-unlist(m$files,use.names=TRUE);stopifnot(all(unname(tools::md5sum(file.path(p,names(ledger))))==unname(ledger)))
 files_n[[site]]<-length(ledger)
 a<-jsonlite::fromJSON(file.path(p,'federated/apache_adapted/audit.json'));stopifnot(a$status=='passed')
}
m<-jsonlite::fromJSON(file.path(root,'reports/site_supplement/collection_manifest.json'));stopifnot(m$apache_adapted_extension_complete,m$patient_level_inputs_read==FALSE)
# Protect the original primary state/output definitions by their prior ledger.
l<-jsonlite::fromJSON(file.path(ap_primary(),'reports/frozen_artifact_ledger.json'),simplifyVector=FALSE)$ledger
for(z in l)stopifnot(unname(tools::md5sum(z$path))==z$md5)
sources<-c(list.files('diagnostics/apache_adapted',full.names=TRUE),list.files('diagnostics/site_return',full.names=TRUE),list.files('code/lib',full.names=TRUE,pattern='[.]R$'),list.files('diagnostics/sofa_completeness',full.names=TRUE,pattern='[.]R$'),'code/02_derive.R','code/04_supplement.R','diagnostics/revision20261022/run.sh','tests/apache_adapted.R','tests/site_supplement.R','tests/site_return_collection.R')
dest<-file.path(root,'provenance/integration_source');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
for(f in sources){to<-file.path(dest,f);dir.create(dirname(to),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(f,to,overwrite=FALSE))}
fed_json(list(status='passed',sites=ap_sites(),aggregate_files_by_site=files_n,cohort_n=sum(counts$n[counts$cohort=='all_case_missing_indicators']),adapted_complete=sum(counts$n[counts$cohort=='complete_hb_proxy']),matched_SOFA=sum(counts$n[counts$cohort=='matched_complete_SOFA']),APACHE_proxy_labeled=TRUE,not_standard_II_or_IV=TRUE,site_return_hashes_valid=TRUE,collector_spec_and_consensus_checks_passed=TRUE,original_primary_artifacts_unchanged=TRUE,source_md5=as.list(tools::md5sum(sources))),file.path(root,'reports/final_integration_audit.json'))
cat('APACHE extension and site-return integration audits passed\n')
