source('diagnostics/revision20261022/lib.R')
root<-'output/run_three_site_refresh_20261026'
# Preserve the prior scientific state, including the old revision analyses.
ledger<-jsonlite::fromJSON('output/run_federated_rehearsal_20261022/reports/frozen_artifact_ledger.json',simplifyVector=FALSE)$ledger
for(v in ledger)stopifnot(unname(tools::md5sum(v$path))==v$md5)
for(site in c('UCMC','NU','RUSH')){
 q<-jsonlite::fromJSON(file.path(old_root(),'sites',site,'export/quality/audit.json'));r<-jsonlite::fromJSON(file.path(rev_root(),'sites',site,'export/revision_audit.json'));c<-jsonlite::fromJSON(file.path(root,'sites',site,'export/comparison_audit.json'))
 stopifnot(q$status=='passed',r$status=='passed',c$status=='passed',c$same_landmark_members)
 if(site!='UCMC')stopifnot(c$source_panel_unchanged)
}
for(v in c('full','loo_UCMC','loo_NU','loo_RUSH')){
 fresh<-jsonlite::fromJSON(file.path(old_root(),'coordinator',v,'frozen_manifest.json'));previous<-jsonlite::fromJSON(file.path('output/run_federated_rehearsal_20261022/coordinator',v,'frozen_manifest.json'))
 stopifnot(fresh$consensus_id!=previous$consensus_id)
}
revaudit<-jsonlite::fromJSON(file.path(rev_root(),'reports/revision_audit.json'));stopifnot(revaudit$status=='passed',revaudit$frozen_consensus_versions==32,revaudit$local_outcome_passes==96)
for(site in c('UCMC','NU','RUSH'))stopifnot(jsonlite::fromJSON(file.path(rev_root(),'sensitivities/absorbing_v2/sites',site,'export/absorbing_audit.json'))$status=='passed')
files<-c(list.files('diagnostics/refresh20261026',recursive=TRUE,full.names=TRUE),list.files('diagnostics/federated',recursive=TRUE,full.names=TRUE),list.files('diagnostics/revision20261022',recursive=TRUE,full.names=TRUE),list.files('code/lib',pattern='[.]R$',full.names=TRUE),'tests/refresh_paths.R')
for(f in files){to<-file.path(root,'provenance/source',f);dir.create(dirname(to),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(f,to,overwrite=FALSE))}
writeLines(capture.output(sessionInfo()),file.path(root,'provenance/sessionInfo.txt'))
fed_json(list(status='passed',date='2026-10-26',active_sites=c('UCMC','NU','RUSH'),MIMIC_included=FALSE,new_UCMC_full_source_run=TRUE,NU_RUSH_source_panels_unchanged=TRUE,all_3_sites_discovery_and_outcome_models_refit=TRUE,new_primary_consensus_versions=4,sensitivity_consensus_versions=32,primary_and_sensitivity_outcome_passes=108,absorbing_descriptive_only=TRUE,all_landmark_members_retained=TRUE,original_primary_frozen_artifacts_unchanged=TRUE,reviewer_audit=revaudit$status,source_md5=as.list(tools::md5sum(files)),privacy='Logical federation rehearsal; linked aggregate package requires local disclosure review'),file.path(root,'reports/final_audit.json'))
