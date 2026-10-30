source('diagnostics/pooled_survival/lib.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)<2)stop('Usage: central_grid.R NEW_OUTPUT SITE_TIME_SHARE [...]')
out<-a[1];shares<-a[-1];if(dir.exists(out))stop('Choose fresh central grid output')
meta<-lapply(shares,function(p)jsonlite::fromJSON(file.path(p,'manifest.json')))
sites<-vapply(meta,function(m)m$site,character(1));stopifnot(!anyDuplicated(sites),length(unique(vapply(meta,function(m)m$consensus_id,character(1))))==1)
for(i in seq_along(shares))for(k in names(meta[[i]]$files))stopifnot(unname(tools::md5sum(file.path(shares[i],k)))==meta[[i]]$files[[k]])
events<-dplyr::bind_rows(lapply(shares,function(p)readr::read_csv(file.path(p,'event_duration_counts.csv'),show_col_types=FALSE)))
cohorts<-dplyr::bind_rows(lapply(shares,function(p)readr::read_csv(file.path(p,'cohort_counts.csv'),show_col_types=FALSE)))
keys<-unique(cohorts[c('curve','endpoint')]);grids<-list()
for(i in seq_len(nrow(keys))){
 k<-keys[i,];t<-events$time_tick[events$curve==k$curve]
 ticks<-sort(unique(c(ps_tick(c(0,7,14,21,28)),t)))
 grids[[i]]<-tibble::tibble(curve=k$curve,endpoint=k$endpoint,time_tick=ticks,time=ps_time(ticks))
}
grid<-dplyr::bind_rows(grids);grid$grid_id<-seq_len(nrow(grid))
ps_export(grid,file.path(out,'shared_grid.csv'));ps_export(cohorts,file.path(out,'site_cohort_counts.csv'))
jsonlite::write_json(list(status='passed',sites=sites,consensus_id=meta[[1]]$consensus_id,grid_md5=unname(tools::md5sum(file.path(out,'shared_grid.csv'))),patient_level_inputs_read=FALSE),file.path(out,'manifest.json'),auto_unbox=TRUE,pretty=TRUE)
cat('Shared event-time grid constructed from aggregate registries\n')
