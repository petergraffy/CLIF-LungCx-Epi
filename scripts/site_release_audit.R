# Read-only source/schema reconciliation on the ORIGINAL Stage 1 checkout.
# Does not load patient rows, rerun analyses, or modify the Stage 1 return.
source('scripts/runtime.R')
a<-commandArgs(trailingOnly=TRUE)
if(length(a)<2L)stop('Usage: Rscript --vanilla /path/to/site_release_audit.R CONFIG REFERENCE_JSON [OUTPUT]')
if(!requireNamespace('jsonlite',quietly=TRUE)||!requireNamespace('digest',quietly=TRUE))stop('Use the project R library containing jsonlite and digest.')
readj<-function(p)jsonlite::fromJSON(p,simplifyVector=FALSE)
writej<-function(x,p)jsonlite::write_json(x,p,auto_unbox=TRUE,pretty=TRUE,null='null',digits=16)
cfg<-readj(a[1]);ref<-readj(a[2]);repo<-normalizePath(path.expand(cfg$repo),winslash='/',mustWork=TRUE)
if(repo!=normalizePath(getwd(),winslash='/'))stop('Run from the original site repository and use its original Stage 1 config.')
site<-cfg$site_name;run<-cfg$analysis$run_id
if(is.null(site)||is.null(run)||!grepl('^[A-Za-z0-9_-]+$',site)||!grepl('^[A-Za-z0-9_-]+$',run))stop('Expected original site_name and analysis.run_id.')
return<-file.path(repo,'output',paste0('run_',site,'_',run,'_federated'),'stage1_return')
m<-readj(file.path(return,'return_manifest.json'))
if(!identical(m$site,site)||m$stage!='stage1')stop('Original Stage 1 return/site mismatch.')
output<-if(length(a)>=3L)path.expand(a[3]) else file.path(repo,'output',paste0('run_',site,'_',run,'_release_audit'))
if(file.exists(output))stop('Audit output already exists; preserve it and choose a fresh output directory.')
dir.create(output,recursive=TRUE);output<-normalizePath(output,winslash='/')
files<-system2(Sys.which('git'),c('ls-files','--','code','diagnostics','scripts','renv.lock','config/waterfall_requirements.txt'),stdout=TRUE)
if(!is.null(attr(files,'status')))stop('Git tracked source inventory unavailable.')
files<-sort(files[grepl('^(code|diagnostics|scripts)/.*[.](R|py|sh|bat)$',files)|files%in%c('renv.lock','config/waterfall_requirements.txt')],method='radix')
if(any(!file.exists(files)))stop('Tracked source is missing; restore the original checkout before auditing.')
textmd5<-function(path){b<-readBin(path,'raw',n=file.info(path)$size);digest::digest(charToRaw(gsub('\r\n','\n',rawToChar(b),fixed=TRUE)),algo='md5',serialize=FALSE)}
source<-as.list(setNames(vapply(files,textmd5,character(1)),files))
reference<-ref$source_files_md5
if(is.null(reference))stop('Reference JSON must contain source_files_md5.')
changed<-sort(union(setdiff(names(source),names(reference)),names(reference)[vapply(names(reference),function(k)!identical(source[[k]],reference[[k]]),logical(1))]),method='radix')
# Only tracked project implementation files; configs, tables and RDS caches are never copied.
for(f in intersect(changed,files)){
 to<-file.path(output,'source_review',f);dir.create(dirname(to),recursive=TRUE,showWarnings=FALSE)
 if(!file.copy(f,to))stop('Unable to copy changed project source: ',f)
}
packet<-lapply(names(m$files),function(rel){
 if(grepl('(^/|(^|/)\\.\\.(/|$))',rel))stop('Unsafe return inventory path.')
 p<-file.path(return,rel);exists<-file.exists(p);ok<-exists&&identical(unname(tools::md5sum(p)),m$files[[rel]])
 list(file=rel,status=if(!exists)'missing' else if(!ok)'hash_mismatch' else 'verified')
})
package_status<-lapply(names(ref$r_packages),function(p){
 found<-tryCatch(as.character(utils::packageVersion(p)),error=function(e)NA_character_)
 list(package=p,expected=ref$r_packages[[p]],observed=found,match=identical(found,ref$r_packages[[p]]))
})
path<-normalizePath(path.expand(cfg$tables_path),winslash='/',mustWork=TRUE)
tables<-list.files(path,recursive=TRUE,full.names=TRUE,pattern='[.]parquet$',ignore.case=TRUE)
bases<-tolower(sub('^clif_','',tools::file_path_sans_ext(basename(tables)),ignore.case=TRUE))
schema<-list();errors<-list()
contract<-m$timestamp_contract
if(is.null(contract))contract<-cfg$timestamp_contract
if(is.null(contract))contract<-cfg$timestamps
for(table in c('patient','hospitalization','adt','respiratory_support','vitals','labs','medication_admin_continuous','patient_assessments','code_status')){
 hit<-tables[bases==table]
 if(length(hit)!=1L){errors[[table]]<-paste('expected one parquet file; found',length(hit));next}
 fields<-tryCatch(arrow::open_dataset(hit)$schema$fields,error=function(e)e)
 if(inherits(fields,'error')){errors[[table]]<-conditionMessage(fields);next}
 for(field in fields){
  if(!grepl('_dttm$',field$name,ignore.case=TRUE))next
  type<-field$type;tz<-if(inherits(type,c('Timestamp','TimestampType')))type$timezone() else ''
  declared<-if(is.null(contract$source_timezone))'' else contract$source_timezone
  schema[[length(schema)+1L]]<-data.frame(table=table,column=field$name,arrow_type=type$ToString(),timezone_aware=nzchar(tz),schema_timezone=tz,declared_source_timezone=declared,review=if(nzchar(tz))'timezone_encoded_in_schema' else if(nzchar(declared))'naive_with_declared_timezone_review_parser' else 'naive_timezone_requires_site_confirmation')
 }
}
if(length(schema))utils::write.csv(do.call(rbind,schema),file.path(output,'timestamp_schema.csv'),row.names=FALSE,na='')
summary<-list(status='audit_complete_requires_coordinator_review',site=site,run_id=run,original_stage1_manifest_md5=unname(tools::md5sum(file.path(return,'return_manifest.json'))),original_stage1_source_hash=m$source_hash,original_timestamp_contract=m$timestamp_contract,current_canonical_source_hash=digest::digest(source,algo='sha256'),current_source_files_md5=source,source_files_differing_from_reference=changed,git_commit=system2(Sys.which('git'),c('rev-parse','HEAD'),stdout=TRUE),r_version=as.character(getRversion()),expected_r_version=ref$r_version,package_status=package_status,timestamp_schema_errors=errors,packet_integrity=packet,source_table_rows_or_patient_rows_read=FALSE,source_schema_metadata_read=TRUE,original_return_or_caches_modified=FALSE,release_approved=FALSE)
writej(summary,file.path(output,'audit_summary.json'))
exports<-list.files(output,recursive=TRUE,full.names=TRUE)
writej(list(stage='source_timestamp_audit',site=site,files=as.list(setNames(unname(tools::md5sum(exports)),substring(exports,nchar(output)+2)))),file.path(output,'audit_manifest.json'))
cat('Return this separate audit directory; preserve the original Stage 1 packet and caches:\n',output,'\nSource-code differences, schemas and package versions only. No patient rows exported.\n',sep='')
