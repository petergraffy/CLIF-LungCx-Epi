# Coordinator: reads ONLY returned aggregate folders; never site-local caches.
source('code/lib/load.R')
args<-commandArgs(trailingOnly=TRUE);if(length(args)<2)stop('Usage: collect.R NEW_OUTPUT SHARE_DIR [SHARE_DIR ...]')
out<-args[1];returns<-args[-1];if(dir.exists(out))stop('Choose a fresh collection output')
meta<-lapply(returns,function(p)jsonlite::fromJSON(file.path(p,'site_return_manifest.json')))
sites<-vapply(meta,function(m)m$site,character(1));stopifnot(!anyDuplicated(sites))
# Validate declared returned files, never follow arbitrary source paths in a manifest.
for(i in seq_along(returns)){
 ledger<-unlist(meta[[i]]$files,use.names=TRUE)
 stopifnot(length(ledger)>0)
 relative<-sub('^.*/share/','',names(ledger))
 stopifnot(!any(grepl('(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$))',relative)))
 actual<-list.files(returns[i],recursive=TRUE,full.names=FALSE)
 actual<-setdiff(actual,'site_return_manifest.json')
 stopifnot(setequal(actual,relative),all(grepl('[.](csv|json)$',actual)))
 observed<-unname(tools::md5sum(file.path(returns[i],relative)))
 stopifnot(identical(observed,unname(as.character(ledger))))
 required_model<-c('federated/primary/discovery/silhouette.csv',
 'federated/primary/outcomes/full/model_metadata.csv','federated/scores/transport_metrics.csv',
 'federated/absorbing/early_events.csv','federated/model_PH_and_denominators.csv','federated/model_PH_terms.csv','federated/exploratory_piecewise_associations.csv','federated/exploratory_piecewise_metadata.csv',paste0('federated/sensitivities/',
 c('h24','h48','onset','support_only','low_flow','carry0','carry12','raw_support'),'/outcomes/full/model_metadata.csv'))
 stopifnot(all(required_model%in%actual))
 for(f in actual[grepl('[.]csv$',actual)]){
  header<-names(readr::read_csv(file.path(returns[i],f),n_max=0,show_col_types=FALSE))
  stopifnot(!any(header%in%c('patient_id','hospitalization_id','t0','event_time','recorded_dttm','admission_dttm','discharge_dttm','sequence')))
 }

}

consensus<-vapply(meta,function(m)m$primary_consensus_id,character(1));stopifnot(length(unique(consensus))==1)
dir.create(out,recursive=TRUE)
files<-lapply(returns,function(p)list.files(file.path(p,'supplement'),pattern='[.]csv$',full.names=FALSE))
required<-c('return_file_inventory.csv','sofa_original_components.csv','sofa_recovery_components.csv','sofa_source_audit_status.csv','site_reporting_metadata.csv')
stopifnot(all(vapply(files,function(f)all(required%in%f),logical(1))))
for(f in sort(unique(unlist(files)))){
 rows<-list()
 empty_template<-NULL
 for(i in seq_along(returns)){
  path<-file.path(returns[i],'supplement',f)
  if(!file.exists(path))next
  x<-readr::read_csv(path,show_col_types=FALSE)
  stopifnot(!any(names(x)%in%c('patient_id','hospitalization_id','t0','event_time','recorded_dttm','admission_dttm','discharge_dttm','sequence')))
  if(nrow(x))stopifnot(all(x$site==sites[i]))
  if(nrow(x))rows[[sites[i]]]<-x else empty_template<-x
 }
 # Keep site rows separate; suppressed cells are not imputed or pooled.
 combined<-if(length(rows))dplyr::bind_rows(rows) else empty_template
 readr::write_csv(combined,file.path(out,f),na='')
}
for(f in c('model_PH_and_denominators.csv','model_PH_terms.csv','exploratory_piecewise_associations.csv','exploratory_piecewise_metadata.csv')){
 rows<-lapply(seq_along(returns),function(i){z<-readr::read_csv(file.path(returns[i],'federated',f),show_col_types=FALSE);if(nrow(z))stopifnot(all(z$site==sites[i]));z})
 readr::write_csv(dplyr::bind_rows(rows),file.path(out,f),na='')
}
apache_flags<-vapply(meta,function(m)isTRUE(m$apache_adapted_extension_complete),logical(1))
if(any(apache_flags)){
 if(!all(apache_flags))stop('Mixed adapted-APACHE return versions; collect consistent site phases')
 component_manifests<-lapply(returns,function(p)jsonlite::fromJSON(file.path(p,'federated/apache_adapted/component_manifest.json')))
 stopifnot(length(unique(vapply(component_manifests,function(m)m$spec_hash,character(1))))==1)
 stopifnot(all(vapply(component_manifests,function(m)m$primary_consensus_id,character(1))==consensus[1]))
 for(p in returns)stopifnot(jsonlite::fromJSON(file.path(p,'federated/apache_adapted/audit.json'))$status=='passed')
 dir.create(file.path(out,'apache_adapted'),recursive=TRUE)
 for(f in c('component_completeness.csv','component_counts.csv','cohort_counts.csv','provenance_counts.csv','missing_patterns.csv','source_metadata_quality.csv','model_metadata.csv','transport_metrics.csv','paired_increments.csv','decision_curves.csv','calibration_bins.csv','evaluation_metadata.csv')){
  rows<-lapply(seq_along(returns),function(i){z<-readr::read_csv(file.path(returns[i],'federated/apache_adapted',f),show_col_types=FALSE);if(nrow(z))stopifnot(all(z$site==sites[i]));z})
  readr::write_csv(dplyr::bind_rows(rows),file.path(out,'apache_adapted',f),na='')
 }
}
descriptor_flags<-vapply(meta,function(m)isTRUE(m$publication_descriptives_complete),logical(1))
if(any(descriptor_flags)){
 if(!all(descriptor_flags))stop('Mixed publication-descriptor return phases')
 for(i in seq_along(returns)){
  dm<-jsonlite::fromJSON(file.path(returns[i],'federated/descriptives/manifest.json'))
  stopifnot(dm$site==sites[i],dm$consensus_id==consensus[1],dm$status=='passed')
 }
 rc<-system2(file.path(R.home('bin'),'Rscript'),c('diagnostics/descriptives/collect.R',shQuote(file.path(out,'publication_tables')),shQuote(returns)))
 if(rc!=0)stop('Publication-table collection failed')
}
curve_flags<-vapply(meta,function(m)isTRUE(m$daily_competing_event_curves_complete),logical(1))
if(any(curve_flags)){
 if(!all(curve_flags))stop('Mixed daily competing-event curve return phases')
 rows<-lapply(seq_along(returns),function(i){
  x<-readr::read_csv(file.path(returns[i],'federated/primary/outcomes/full/cumulative_incidence_daily.csv'),show_col_types=FALSE)
  stopifnot(all(x$site==sites[i]),all(x$consensus_id==consensus[1]))
  x
 })
 readr::write_csv(dplyr::bind_rows(rows),file.path(out,'cumulative_incidence_daily.csv'),na='')
}
survival_flags<-vapply(meta,function(m)isTRUE(m$survival_tdauc_complete),logical(1))
if(any(survival_flags)){
 if(!all(survival_flags))stop('Mixed survival/tdAUC return phases')
 for(i in seq_along(returns)){
  sm<-jsonlite::fromJSON(file.path(returns[i],'federated/survival_tdauc/manifest.json'))
  stopifnot(sm$site==sites[i],sm$primary_consensus_id==consensus[1],sm$status=='passed')
 }
 dir.create(file.path(out,'survival_tdauc'))
 for(f in c('survival_curves.csv','survival_denominators.csv','tdauc.csv','tdroc.csv')){
  rows<-lapply(seq_along(returns),function(i){x<-readr::read_csv(file.path(returns[i],'federated/survival_tdauc',f),show_col_types=FALSE);stopifnot(all(x$site==sites[i]));x})
  readr::write_csv(dplyr::bind_rows(rows),file.path(out,'survival_tdauc',f),na='')
 }
}
pool_flags<-vapply(meta,function(m)isTRUE(m$pooled_survival_discharge_complete),logical(1))
if(identical(Sys.getenv('CLIF_REQUIRE_POOLED_SURVIVAL'),'true')&&!all(pool_flags))stop('Final publication return requires pooled survival/discharge exports from every site')
if(any(pool_flags)){
 if(!all(pool_flags))stop('Mixed pooled survival/discharge return phases')
 grid_hashes<-character(length(returns))
 for(i in seq_along(returns)){
  base<-file.path(returns[i],'federated/pooled_survival')
  tm<-jsonlite::fromJSON(file.path(base,'times/manifest.json'))
  rm<-jsonlite::fromJSON(file.path(base,'risksets/manifest.json'))
  stopifnot(tm$site==sites[i],rm$site==sites[i],tm$consensus_id==consensus[1],rm$consensus_id==consensus[1],tm$status=='passed',rm$status=='passed',rm$bootstrap==200,rm$local_AJ_reconstruction_passed,rm$local_KM_reconstruction_passed)
  for(phase in c('times','risksets')){
   p<-file.path(base,phase);em<-jsonlite::fromJSON(file.path(p,'manifest.json'))
   for(k in names(em$files))stopifnot(unname(tools::md5sum(file.path(p,k)))==em$files[[k]])
  }
  q<-readr::read_csv(file.path(base,'times/cohort_counts.csv'),show_col_types=FALSE)
  stopifnot(all(c('death','death_or_hospice','live_discharge','any_exit')%in%q$endpoint))
  grid_hashes[i]<-rm$grid_md5
 }
 stopifnot(length(unique(grid_hashes))==1)
}
cox_flags<-vapply(meta,function(m)isTRUE(m$four_endpoint_cox_complete),logical(1))
if(identical(Sys.getenv('CLIF_REQUIRE_FOUR_ENDPOINT_COX'),'true')&&!all(cox_flags))stop('Final return requires four-endpoint Cox exports from every site')
if(any(cox_flags)){
 if(!all(cox_flags))stop('Mixed four-endpoint Cox return phases')
 extensions<-file.path(returns,'federated/four_endpoint_cox')
 for(i in seq_along(extensions)){
  cm<-jsonlite::fromJSON(file.path(extensions[i],'manifest.json'))
  stopifnot(cm$site==sites[i],cm$primary_consensus_id==consensus[1],cm$status=='passed',setequal(cm$endpoints,cox_endpoints()))
 }
 rc<-system2(file.path(R.home('bin'),'Rscript'),c('diagnostics/four_endpoint_cox/collect.R',shQuote(file.path(out,'four_endpoint_cox')),shQuote(extensions)))
 if(rc!=0)stop('Four-endpoint Cox collection failed')
}
coverage<-dplyr::bind_rows(lapply(seq_along(returns),function(i)tibble::tibble(site=sites[i],family=unique(unlist(files)),status=ifelse(unique(unlist(files))%in%files[[i]],'received','not_returned'))))
readr::write_csv(coverage,file.path(out,'collection_coverage.csv'))
jsonlite::write_json(list(status='aggregate_collection_complete',sites=sites,primary_consensus_id=consensus[1],site_count=length(sites),counts_pooled=FALSE,suppressed_values_preserved=TRUE,patient_level_inputs_read=FALSE,author_confirmation_required=TRUE,apache_adapted_extension_complete=all(apache_flags),publication_descriptives_complete=all(descriptor_flags),daily_competing_event_curves_complete=all(curve_flags),survival_tdauc_complete=all(survival_flags),pooled_survival_discharge_complete=all(pool_flags),four_endpoint_cox_complete=all(cox_flags)),file.path(out,'collection_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
