# First local pass: aggregate event-duration registry, with no calendar dates/IDs.
source('diagnostics/pooled_survival/lib.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=3)stop('Usage: site_times.R SITE PRIMARY_ROOT NEW_OUTPUT')
site<-a[1];primary<-a[2];out<-a[3];if(dir.exists(out))stop('Choose fresh local grid output')
x<-ps_data(site,primary);rows<-list();meta<-list()
for(key in names(x$sets))for(endpoint in c('death','death_or_hospice','live_discharge','any_exit')){
 d<-x$sets[[key]];event<-ps_event(d,endpoint)
 z<-tibble::tibble(time_tick=ps_tick(d$time[event==1&d$time<=28])) |> dplyr::count(time_tick,name='events') |> dplyr::mutate(time=ps_time(time_tick))
 rows[[paste(key,endpoint)]]<-z |> dplyr::mutate(curve=key,endpoint=endpoint,site=site)
 meta[[paste(key,endpoint)]]<-tibble::tibble(site=site,curve=key,endpoint=endpoint,n=nrow(d),max_followup=if(nrow(d))max(d$time) else 0,unknown_outcome_n=sum(!d$valid_outcome))
}
ps_export(dplyr::bind_rows(rows),file.path(out,'event_duration_counts.csv'))
ps_export(dplyr::bind_rows(meta),file.path(out,'cohort_counts.csv'))
files<-list.files(out,full.names=TRUE)
jsonlite::write_json(list(status='passed',site=site,consensus_id=jsonlite::fromJSON(file.path(primary,'coordinator/full/frozen_manifest.json'))$consensus_id,
 durations='Aggregate relative follow-up durations in days; no absolute timestamps or identifiers',files=as.list(setNames(unname(tools::md5sum(files)),basename(files)))),file.path(out,'manifest.json'),auto_unbox=TRUE,pretty=TRUE)
cat(site,'aggregate event-time registry exported\n')
