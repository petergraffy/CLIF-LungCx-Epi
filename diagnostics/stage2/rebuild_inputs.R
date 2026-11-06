source('scripts/deployment_lib.R');source('diagnostics/stage2/rebuild_lib.R')
a<-commandArgs(trailingOnly=TRUE);cfg<-deploy_config(Sys.getenv('CLIF_RUN_CONFIG'))
if(length(a)!=1L||!a%in%c('raw','waterfall'))stop('Usage: rebuild_inputs.R raw|waterfall')
if(a=='raw'){
 rows<-list();files<-list.files(cfg$tables_path,pattern='[.]parquet$',recursive=TRUE,full.names=TRUE,ignore.case=TRUE)
 for(f in files){schema<-arrow::open_dataset(f)$schema
  for(ts in grep('_dttm$',schema$names,value=TRUE)){
   type<-schema$GetFieldByName(ts)$type$ToString();rows[[paste(f,ts)]]<-tibble::tibble(table=tolower(sub('^clif_','',tools::file_path_sans_ext(basename(f)))),column=ts,arrow_type=type,handling=if(grepl('tz=',type,fixed=TRUE))'preserve_absolute_instant' else paste0('interpret_naive_in_',cfg$source_timezone),dst_policy='ambiguous_or_nonexistent_local_time_is_missing')
  }
 }
 dir<-file.path(dirname(cfg$work_root),'input_audit');dir.create(dir,recursive=TRUE,showWarnings=FALSE);readr::write_csv(dplyr::bind_rows(rows),file.path(dir,'source_timestamp_schema.csv'))
 s2_raw_inputs(cfg)
}else s2_waterfall_inputs(cfg,Sys.getenv('CLIF_WATERFALL_OUTPUT'))
