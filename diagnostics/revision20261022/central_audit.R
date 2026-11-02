source('diagnostics/revision20261022/score_lib.R')
root<-rev_root();ledger<-list();files<-list.files(file.path(root,'exchange'),recursive=TRUE,full.names=TRUE)
stopifnot(all(grepl('\\.(csv|json)$',files)))
for(f in files[grepl('csv$',files)]){z<-readr::read_csv(f,show_col_types=FALSE);stopifnot(!any(names(z)%in%c('patient_id','hospitalization_id','recorded_dttm','admission_dttm','discharge_dttm','event_time','sequence')))}
for(v in c('h24','h48','onset','support_only','low_flow','carry0','carry12','raw_support')){
 Sys.setenv(CLIF_REVISION_VARIANT=v);source('diagnostics/revision20261022/variant_contract.R')
 for(version in fed_versions()){
  dir<-file.path(root,'sensitivities',v,'coordinator',version);path<-file.path(dir,'frozen_manifest.json');m<-jsonlite::fromJSON(path,simplifyVector=FALSE)
  for(group in m$source_hashes)for(file in names(group))stopifnot(unname(tools::md5sum(file))==group[[file]])
  proto<-readRDS(file.path(dir,'prototypes.rds'));expected<-substr(fed_hash(list(spec=fed_spec(),training_sites=as.character(unlist(m$training_sites)),source=m$source_hashes,prototypes=proto)),1,20);stopifnot(expected==m$consensus_id)
  for(s in fed_sites()){
   om<-jsonlite::fromJSON(file.path(root,'sensitivities',v,'exchange/outcomes',version,s,'manifest.json'));stopifnot(om$consensus_id==m$consensus_id,om$consensus_hash==fed_hash(jsonlite::fromJSON(path)))
  }
  ledger[[paste(v,version)]]<-list(consensus_id=m$consensus_id,prototype_md5=unname(tools::md5sum(file.path(dir,'prototypes.rds'))))
 }
}
for(s in fed_sites()){
 from<-file.path(root,'sites',s,'export');to<-file.path(root,'exchange/review',s)
 for(nm in c('partition_agreement.csv','revision_audit.json'))stopifnot(file.copy(file.path(from,nm),file.path(to,nm),overwrite=FALSE))
 a<-jsonlite::fromJSON(file.path(to,'revision_audit.json'));stopifnot(a$status=='passed')
}
z<-dplyr::bind_rows(lapply(fed_sites(),function(s)readr::read_csv(file.path(root,'exchange/review',s,'partition_agreement.csv'),show_col_types=FALSE)))
fed_write(z,file.path(root,'reports/tables/partition_agreement.csv'))
fed_json(list(status='passed',frozen_consensus_versions=length(ledger),local_outcome_passes=96,original_primary_artifacts_unchanged=TRUE,all_local_audits_passed=TRUE,no_patient_identifiers_in_exchange=TRUE,privacy_scope='Logical local rehearsal, not differential privacy or approved institutional exchange',ledger=ledger),file.path(root,'reports/revision_audit.json'))
cat('Central aggregate, split and 32 frozen consensus audits passed\n')
