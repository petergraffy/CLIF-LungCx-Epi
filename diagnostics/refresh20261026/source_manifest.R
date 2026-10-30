source('code/lib/load.R')
cfg<-read_run_config('output/run_three_site_refresh_20261026/local/ucmc_config.json');base<-file.path(cfg$repo,paste0('output/run_UCMC_',cfg$analysis$run_id));p<-readRDS(file.path(base,'local/prepared.rds'))
files<-file.path(cfg$tables_path,paste0('clif_',names(p$tables),'.parquet'));files<-files[file.exists(files)]
rows<-lapply(files,function(f)tibble::tibble(table=basename(f),bytes=file.info(f)$size,md5=unname(tools::md5sum(f))))
readr::write_csv(dplyr::bind_rows(rows),'output/run_three_site_refresh_20261026/local/ucmc_source_fingerprints.csv')
jsonlite::write_json(list(site='UCMC',source_path=cfg$tables_path,fresh_baseline_root=base,run_id=cfg$analysis$run_id,non_UCMC='NU/Rush validated Oct22 source caches reused; all discovery/outcome/sensitivity fits redone',source_code=as.list(tools::md5sum(list.files('code',pattern='[.]R$',recursive=TRUE,full.names=TRUE)))), 'output/run_three_site_refresh_20261026/local/input_manifest.json',pretty=TRUE,auto_unbox=TRUE)
