# Run on one site's machine after its frozen primary/reviewer passes.
source('code/lib/load.R')
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=4)stop('Usage: Rscript diagnostics/site_return/capture.R SITE PRIMARY_ROOT REVISION_ROOT NEW_RETURN_ROOT')
site<-args[1];primary_root<-normalizePath(args[2]);revision_root<-normalizePath(args[3]);out<-args[4]
if(dir.exists(out))stop('Choose a fresh site-return root; completed returns are immutable')
home<-file.path(primary_root,'sites',site);review<-file.path(revision_root,'sites',site)
discovery<-readRDS(file.path(home,'local/discovery.rds'));input<-jsonlite::fromJSON(file.path(discovery$input_root,'local/input_manifest.json'))
stopifnot(input$site==site)
baseline<-input$baseline_root;prepared<-readRDS(file.path(baseline,'local/prepared.rds'))
cfg<-list(site_name=site,repo=normalizePath('.'),tables_path=Sys.getenv(paste0('CLIF_SITE_TABLES_',site),unset=dirname(input$source_file)),file_type='parquet',analysis=prepared$parameters)
dir.create(file.path(out,'local'),recursive=TRUE)
config<-file.path(out,'local/site_config.json');jsonlite::write_json(cfg,config,auto_unbox=TRUE,pretty=TRUE)
Sys.setenv(CLIF_SUPPLEMENT_OUTPUT=out)
capture_site_supplement(cfg,baseline,config)
# The baseline supplements never relabel baseline clusters as shared phenotypes.
# Collect shared-phenotype/model tables only from the actual frozen site passes.
shared<-file.path(out,'share/federated');dir.create(shared,recursive=TRUE)
copy_aggregate<-function(src,dest){
 if(!file.exists(src))stop('Required site return missing: ',basename(src))
 ext<-tolower(tools::file_ext(src));stopifnot(ext%in%c('csv','json'))
 if(ext=='csv'){
  x<-readr::read_csv(src,show_col_types=FALSE)
  stopifnot(!any(names(x)%in%c('patient_id','hospitalization_id','t0','event_time','recorded_dttm','admission_dttm','discharge_dttm','sequence')))
 } else {
  x<-jsonlite::fromJSON(src,simplifyVector=FALSE)
  # Source manifests can contain private file paths: share hashes/specification only.
  scrub<-function(z){if(!is.list(z))return(z);for(k in names(z)){
   if(grepl('path|input_files|discovery_inputs|source_file|source_root',k,ignore.case=TRUE))z[[k]]<-NULL else z[[k]]<-scrub(z[[k]])};z}
  x<-scrub(x);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE);jsonlite::write_json(x,dest,auto_unbox=TRUE,pretty=TRUE,null='null');return(invisible(TRUE))
 }
 dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(src,dest))
}
# Explicit allowlists prevent a local RDS/cache/prediction file from entering returns.
primary<-list(discovery=c('silhouette.csv','clusters.csv','fingerprints.csv','manifest.json'),
 risks='risks_28day.csv',transport=c('metrics.csv','coverage.csv','manifest.json','extended_metrics.csv','extended_manifest.json'),quality='audit.json')
for(family in names(primary))for(f in primary[[family]])copy_aggregate(file.path(home,'export',family,f),file.path(shared,'primary',family,f))
versions<-list.dirs(file.path(home,'export/outcomes'),full.names=FALSE,recursive=FALSE)
stopifnot('full'%in%versions,paste0('loo_',site)%in%versions)
for(v in versions)for(f in c('assignment_counts.csv','associations.csv','association_covariances.csv','model_metadata.csv','cumulative_incidence.csv','prediction_coefficients.csv','manifest.json'))copy_aggregate(file.path(home,'export/outcomes',v,f),file.path(shared,'primary/outcomes',v,f))
review_files<-c('arf_associations.csv','arf_joint_tests.csv','arf_type_profiles.csv','clinical_table.csv','code_source_audit.json','code_source_catalogue.csv','code_transitions.csv','cohort_flow.csv','continuous_denominators.csv','early_events.csv','outcomes_by_code.csv','partition_agreement.csv','phenotype_characteristics.csv','review_manifest.json','revision_audit.json','sensitivity_associations.csv','sensitivity_metadata.csv','severity_availability.csv','severity_source_audit.csv')
for(f in review_files)copy_aggregate(file.path(review,'export',f),file.path(shared,'review',f))
for(f in c('transport_metrics.csv','paired_increments.csv','decision_curves.csv','transport_manifest.json','manifest.json','model_metadata.csv','coefficients.csv'))copy_aggregate(file.path(review,'export/scores',f),file.path(shared,'scores',f))
variants<-c('h24','h48','onset','support_only','low_flow','carry0','carry12','raw_support')
for(a in variants){
 h<-file.path(revision_root,'sensitivities',a,'sites',site,'export')
 for(f in c('silhouette.csv','clusters.csv','fingerprints.csv','manifest.json'))copy_aggregate(file.path(h,'discovery',f),file.path(shared,'sensitivities',a,'discovery',f))
 for(v in versions)for(f in c('assignment_counts.csv','associations.csv','association_covariances.csv','model_metadata.csv','cumulative_incidence.csv','prediction_coefficients.csv','manifest.json'))copy_aggregate(file.path(h,'outcomes',v,f),file.path(shared,'sensitivities',a,'outcomes',v,f))
 for(f in c('metrics.csv','coverage.csv','manifest.json'))copy_aggregate(file.path(h,'transport',f),file.path(shared,'sensitivities',a,'transport',f))
}
h<-file.path(revision_root,'sensitivities/absorbing_v2/sites',site,'export')
for(f in c('early_events.csv','absorbing_manifest.json','absorbing_audit.json'))copy_aggregate(file.path(h,f),file.path(shared,'absorbing',f))
for(f in c('silhouette.csv','clusters.csv','fingerprints.csv','manifest.json'))copy_aggregate(file.path(h,'discovery',f),file.path(shared,'absorbing/discovery',f))
# Waterfall QA has exact local diagnostics; apply this site's disclosure threshold.
f<-file.path(discovery$input_root,'aggregate/waterfall_source_conflict_checks.csv')
if(file.exists(f)){
 x<-readr::read_csv(f,show_col_types=FALSE);x<-protect_aggregate(x,intersect(c('rows','hours','hospitalizations'),names(x)),character(),cfg$analysis$min_cell)
 write_share(x,file.path(out,'local/supplement_package'),'waterfall_source_conflict_checks',cfg)
 copy_aggregate(file.path(out,'local/supplement_package/share/waterfall_source_conflict_checks.csv'),file.path(shared,'waterfall_source_conflict_checks.csv'))
}
# Global PH tests and actual model denominators are checked locally, with safe returns.
x<-readRDS(file.path(home,'local/full/outcomes.rds'))
metadata<-readr::read_csv(file.path(home,'export/outcomes/full/model_metadata.csv'),show_col_types=FALSE)
stopifnot(setequal(metadata$endpoint,cox_endpoints()))
diagnostic<-supplement_model_diagnostics(x,metadata)
saveRDS(diagnostic,file.path(out,'local/model_PH_diagnostic_fits.rds'))
ph<-diagnostic$table;readr::write_csv(ph,file.path(out,'local/model_PH_and_denominators_exact.csv'))
ph<-protect_aggregate(ph,c('n','events_n','non_events_n'),character(),cfg$analysis$min_cell)
write_share(ph,file.path(out,'local/supplement_package'),'model_PH_and_denominators',cfg)
copy_aggregate(file.path(out,'local/supplement_package/share/model_PH_and_denominators.csv'),file.path(shared,'model_PH_and_denominators.csv'))
terms<-protect_aggregate(diagnostic$term_tests,c('n','events_n','non_events_n'),'model_key',cfg$analysis$min_cell)
write_share(terms,file.path(out,'local/supplement_package'),'model_PH_terms',cfg)
copy_aggregate(file.path(out,'local/supplement_package/share/model_PH_terms.csv'),file.path(shared,'model_PH_terms.csv'))
piece<-supplement_piecewise_associations(x$data,cfg$analysis$min_cell)
saveRDS(piece,file.path(out,'local/piecewise_diagnostic_fits.rds'))
write_share(piece$associations,file.path(out,'local/supplement_package'),'exploratory_piecewise_associations',cfg)
meta<-protect_aggregate(piece$metadata,c('candidate_n','fitted_n','events_n','non_events_n'),character(),cfg$analysis$min_cell)
write_share(meta,file.path(out,'local/supplement_package'),'exploratory_piecewise_metadata',cfg)
for(f in c('exploratory_piecewise_associations.csv','exploratory_piecewise_metadata.csv'))copy_aggregate(file.path(out,'local/supplement_package/share',f),file.path(shared,f))

apache_root<-Sys.getenv('CLIF_APACHE_OUTPUT')
apache_present<-nzchar(apache_root)
if(apache_present){
 ap<-file.path(apache_root,'sites',site,'export')
 expected_ap<-c('component_completeness.csv','component_counts.csv','cohort_counts.csv','provenance_counts.csv','missing_patterns.csv','source_metadata_quality.csv','coefficients.csv','model_metadata.csv','transport_metrics.csv','paired_increments.csv','decision_curves.csv','calibration_bins.csv','evaluation_metadata.csv','component_manifest.json','training_manifest.json','transport_manifest.json','audit.json')
 for(f in expected_ap)copy_aggregate(file.path(ap,f),file.path(shared,'apache_adapted',f))
}
# Shared publication Table2 is generated only from the frozen full assignments.
primary_data<-readRDS(file.path(home,'local/full/outcomes.rds'))$data
source('diagnostics/manuscript_figures/cif_lib.R')
daily_cif<-released_cif_grid(primary_data,0:28,cfg$analysis$min_cell) |> dplyr::mutate(site=site,version='full',consensus_id=jsonlite::fromJSON(file.path(primary_root,'coordinator/full/frozen_manifest.json'))$consensus_id)
readr::write_csv(daily_cif,file.path(shared,'primary/outcomes/full/cumulative_incidence_daily.csv'),na='')
ap_components<-if(apache_present)readRDS(file.path(apache_root,'sites',site,'local/components.rds')) else readRDS(file.path(out,'local/apache_adapted/components.rds'))
shared_root<-file.path(out,'local/shared_descriptors')
phenotypes<-rownames(readRDS(file.path(primary_root,'coordinator/full/prototypes.rds')))
consensus_id<-jsonlite::fromJSON(file.path(primary_root,'coordinator/full/frozen_manifest.json'))$consensus_id
capture_descriptive_tables(prepared,primary_data,cfg,shared_root,discovery$panel,phenotypes,ap_components,consensus_id)
for(f in list.files(file.path(shared_root,'share/descriptives'),full.names=TRUE))copy_aggregate(f,file.path(shared,'descriptives',basename(f)))
for(f in c('cohort_exclusions.csv','cohort_by_year.csv','cohort_patient_counts.csv','outcome_timestamp_audit.csv','source_availability.csv','analysis_plan.csv','arf_detection_audit.csv','derivation_status.csv'))copy_aggregate(file.path(baseline,'share',f),file.path(shared,'cohort_audit',f))
files<-list.files(file.path(out,'share'),recursive=TRUE,full.names=TRUE)
stopifnot(!any(grepl('[.](rds|parquet|fst)$',files)))
manifest<-list(status='complete_local_capture_pending_site_review',site=site,run_id=cfg$analysis$run_id,primary_consensus_id=jsonlite::fromJSON(file.path(home,'export/outcomes/full/manifest.json'))$consensus_id,
 versions=versions,primary_cause_specific_cox_endpoints=cox_endpoints(),inferential_sensitivities=variants,absorbing_descriptive_only=TRUE,SOFA_recovery_diagnostic_only_models_use_frozen_original_scores=TRUE,apache_adapted_extension_complete=apache_present,apache_not_standard_II_or_IV=TRUE,publication_descriptives_complete=TRUE,publication_table2_grouping='frozen_shared_primary',daily_competing_event_curves_complete=TRUE,
 source_method='constant-cost local OM; interval fingerprints; frozen per-hospitalization template classification',
 files=as.list(setNames(unname(tools::md5sum(files)),substring(files,nchar(file.path(out,'share'))+2))),capture_source_md5=as.list(tools::md5sum(c('diagnostics/site_return/capture.R','code/lib/site_supplement.R'))),privacy='Return share/ only after local linked-package review; private clinical inputs and diagnostics remain local/')
jsonlite::write_json(manifest,file.path(out,'share/site_return_manifest.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
message('Complete site supplement/return captured: ',file.path(out,'share'))
