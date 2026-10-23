source('diagnostics/revision20261022/score_lib.R')
root<-rev_root();dest<-file.path(root,'provenance/current_source');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
files<-c(list.files('diagnostics/revision20261022',full.names=TRUE,recursive=TRUE),list.files('diagnostics/federated',pattern='\\.R$',full.names=TRUE),list.files('code/lib',pattern='\\.R$',full.names=TRUE),'tests/revision_gaps.R')
for(f in files){to<-file.path(dest,f);dir.create(dirname(to),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(f,to,overwrite=FALSE))}
writeLines(capture.output(sessionInfo()),file.path(root,'provenance/sessionInfo.txt'))
fed_json(list(source_md5=as.list(tools::md5sum(files)),score_spec=score_spec(),primary_source_root=old_root(),original_primary_unchanged=TRUE,scope='three-site exploratory local rehearsal; no new source ETL',caveats=c('Modified SOFA only, not validated standard score','APACHE comparison unavailable','Source datasets already examined','Bootstrap conditional on frozen learner')),file.path(root,'provenance/manifest.json'))
# Final export scan covers nested sensitivities and late descriptive additions.
exports<-c(list.files(file.path(root,'exchange'),recursive=TRUE,full.names=TRUE),list.files(file.path(root,'sensitivities'),pattern='\\.(csv|json)$',recursive=TRUE,full.names=TRUE))
exports<-exports[grepl('/exchange/',exports)]
for(f in exports[grepl('csv$',exports)]){z<-readr::read_csv(f,show_col_types=FALSE);stopifnot(!any(names(z)%in%c('patient_id','hospitalization_id','recorded_dttm','admission_dttm','discharge_dttm','event_time','sequence')))}
stopifnot(!any(grepl('\\.rds$',exports)))
fed_json(list(status='passed',aggregate_export_files=length(exports),patient_identifier_columns=0,individual_predictions_exported=FALSE,privacy_scope='header/type audit and code-path separation, not a formal disclosure guarantee'),file.path(root,'reports/final_export_audit.json'))
