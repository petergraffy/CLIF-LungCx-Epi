# Portable site context. Clinical source paths and caches are never shared.
sofa_context <- function(site) {
 config <- Sys.getenv('CLIF_SOFA_CONFIG')
 if (!nzchar(config)) stop('Set CLIF_SOFA_CONFIG to this site config; no site-specific paths are inferred')
 cfg <- read_run_config(config)
 if (!identical(site,cfg$site_name)) stop('SOFA site/config mismatch')
 baseline <- Sys.getenv('CLIF_SOFA_BASELINE',unset=file.path(cfg$repo,'output',paste0('run_',site,'_',cfg$analysis$run_id)))
 output <- Sys.getenv('CLIF_SOFA_OUTPUT',unset=file.path(baseline,'local','supplement_sofa'))
 for(k in c('local','reports'))dir.create(file.path(output,k),recursive=TRUE,showWarnings=FALSE)
 list(cfg=cfg,roots=setNames(baseline,site),paths=setNames(cfg$tables_path,site),output=output)
}
sofa_source_rows <- function(path, table, ids) {
 files<-list.files(path,full.names=TRUE,recursive=TRUE,pattern='\\.(parquet|csv|fst)$',ignore.case=TRUE)
 bases<-tolower(sub('^clif_','',tools::file_path_sans_ext(basename(files)),ignore.case=TRUE))
 files<-files[bases==table];if(length(files)!=1)stop('Expected one source table: ',table)
 ext<-tolower(tools::file_ext(files))
 if(ext=='parquet'){
  dataset<-arrow::open_dataset(files);x<-dataset |> dplyr::filter(hospitalization_id%in%ids) |> dplyr::collect()
  for(ts in grep('_dttm$',names(x),value=TRUE))x[[ts]]<-source_time(x[[ts]],dataset$schema$GetFieldByName(ts)$type)
 } else {
  x<-if(ext=='csv')readr::read_csv(files,show_col_types=FALSE) else fst::read_fst(files)
  names(x)<-tolower(names(x));x<-x[x$hospitalization_id%in%ids,,drop=FALSE]
 }
 names(x)<-tolower(names(x));x
}
