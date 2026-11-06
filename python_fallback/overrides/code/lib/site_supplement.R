# Site-local supplement capture. Never modifies source tables, scores or fits.
supplement_path <- function(cfg) file.path(cfg$repo,'output',paste0('run_',cfg$site_name,'_',cfg$analysis$run_id))
supplement_write <- function(x, root, name, cfg, counts=character(), partition=character()) {
 dir.create(file.path(root,'local','supplement_exact'),recursive=TRUE,showWarnings=FALSE)
 readr::write_csv(x,file.path(root,'local','supplement_exact',paste0(name,'.csv')),na='')
 if(length(counts))x<-protect_aggregate(x,counts,partition,cfg$analysis$min_cell)
 pkg<-file.path(root,'local','supplement_package');dir.create(file.path(pkg,'share'),recursive=TRUE,showWarnings=FALSE)
 write_share(x,pkg,name,cfg)
}
site_supplement_tables <- function(prepared, analyses, cfg, root=prepared$root) {
 co<-prepared$cohort
 primary<-analyses$h72_primary$data
 stopifnot(!is.null(primary),!anyDuplicated(primary$hospitalization_id))
 for(k in c('local','share'))dir.create(file.path(root,k),recursive=TRUE,showWarnings=FALSE)
 # Bounded windows are documented even when source evidence is absent.
 domains<-c('sofa_resp','sofa_cv','sofa_coag','sofa_liver','sofa_renal','sofa_cns')
 rows<-list();patterns<-list()
 for(scope in c('eligible_cohort','primary_landmark')){
  d<-if(scope=='eligible_cohort')co else primary
  for(v in c(domains,'sofa_total'))rows[[paste(scope,v)]]<-tibble::tibble(cohort=scope,component=v,n=nrow(d),available_n=sum(!is.na(d[[v]])),missing_n=sum(is.na(d[[v]])),method='original_measured_component_modified_SOFA',window='first24h_after_first_ICU_arrival')
  q<-as.data.frame(is.na(d[,domains]));q$pattern<-apply(q,1,function(z)paste(ifelse(z,'missing','observed'),collapse='|'))
  patterns[[scope]]<-tibble::tibble(cohort=scope,pattern=q$pattern) |> dplyr::count(cohort,pattern,name='n')
 }
 supplement_write(dplyr::bind_rows(rows),root,'sofa_original_components',cfg,c('n','available_n','missing_n'),'component')
 supplement_write(dplyr::bind_rows(patterns),root,'sofa_missing_patterns',cfg,'n','cohort')
 # Count documentation, not presumed absence of cancer characteristics/preferences.
 optional<-c('histology','stage','molecular_subtype','ecog','therapy_line','therapy_intent','treatment_type')
 availability<-dplyr::bind_rows(lapply(optional,function(v){present<-v%in%names(primary);n<-if(present)sum(!is.na(primary[[v]])&nzchar(trimws(as.character(primary[[v]])))) else 0L
  tibble::tibble(variable=v,source_present=present,n=nrow(primary),observed_n=n,missing_n=nrow(primary)-n,status=if(!present)'not_supplied' else if(n==0)'supplied_without_values' else 'descriptive_timing_unverified')}))
 supplement_write(availability,root,'oncology_availability',cfg,c('n','observed_n','missing_n'),'variable')
 flow<-co |> dplyr::mutate(reason=dplyr::case_when(!valid_outcome~'unclassifiable_outcome',event_time<=t0+72*3600~paste0('early_',event_type),TRUE~'primary_landmark')) |> dplyr::count(reason,name='n')
 supplement_write(flow,root,'landmark_reconciliation',cfg,'n',character())
 repeat_counts<-tibble::tibble(n=nrow(primary),patients_n=dplyr::n_distinct(primary$patient_id),repeat_hospitalizations_n=sum(!primary$first_hospitalization))
 supplement_write(repeat_counts,root,'repeated_admissions',cfg,names(repeat_counts),character())
 # Every derivation, including an unavailable/failed one, has an explicit status.
 registry<-analysis_specifications(cfg);registry$status<-ifelse(registry$analysis_id%in%names(analyses),'captured','not_available_see_derivation_status')
 supplement_write(registry,root,'analysis_return_inventory',cfg)
 early<-analyses$h72_absorbing$data
 if(!is.null(early)){
  early$early_event<-ifelse(early$event_time<=early$t0+72*3600,early$event_type,'alive_in_hospital_at72h')
  e<-early |> dplyr::count(local_cluster,early_event,name='n') |> dplyr::mutate(group_type=if(is.null(analyses$h72_absorbing$grouping))'local_absorbing_descriptive_only' else 'frozen_absorbing_descriptive_only',analysis_id='h72_absorbing')
  supplement_write(e,root,'early_events_by_absorbing_cluster',cfg,'n','local_cluster')
 }
 # Treatment documentation and terminal outcomes remain separate variables.
 for(v in intersect(c('code_at_admission','code_at_landmark','new_dni','new_comfort','code_changed','any_imv','any_vasopressor','metastatic_cancer','malignant_pleural_effusion','cachexia'),names(primary))){
  q<-primary;q$category<-ifelse(is.na(q[[v]]),'unknown',as.character(q[[v]]))
  z<-q |> dplyr::count(local_cluster,category,name='n') |> dplyr::mutate(variable=v)
  supplement_write(z,root,paste0('clinical_',v),cfg,'n','local_cluster')
 }
 # Annual observation patterns and ARF contribution, computed from actual sequences.
 hourly<-list();arf<-list()
 for(id in names(analyses)){
  a<-analyses[[id]];if(is.null(a$fit$sequence))next
  sq<-a$fit$sequence;d<-a$data[match(a$fit$ids,a$data$hospitalization_id),,drop=FALSE]
  support<-sub('[|].*$','',sq);phys<-sub('^[^|]*[|]','',sq)
  terminal<-support%in%TERMINAL_STATES;dim(terminal)<-dim(sq)
  measures<-list(unknown_support=support=='UNKNOWN',terminal=terminal)
  if(a$fit$field!='support')measures<-c(measures,list(hypox_unknown=grepl('^HU_',phys),hyper_unknown=grepl('_CU$',phys),hypox_positive=grepl('^H1_',phys),hyper_positive=grepl('_C1$',phys)))
  for(v in names(measures)){
   z<-measures[[v]];dim(z)<-dim(sq)
   annual<-tibble::tibble(year=d$admission_year,n=ncol(sq),hours_n=rowSums(z)) |> dplyr::group_by(year) |> dplyr::summarise(panel_hours_n=sum(n),state_hours_n=sum(hours_n),.groups='drop') |> dplyr::mutate(analysis_id=id,measure=v)
   hourly[[paste(id,v)]]<-annual
   if(v%in%c('hypox_positive','hyper_positive'))arf[[paste(id,v)]]<-tibble::tibble(analysis_id=id,local_cluster=d$local_cluster,measure=v,positive=rowSums(z)>0) |> dplyr::count(analysis_id,local_cluster,measure,positive,name='n')
  }
 }
 supplement_write(dplyr::bind_rows(hourly),root,'annual_hourly_evidence',cfg,c('panel_hours_n','state_hours_n'),c('analysis_id','year'))
 supplement_write(dplyr::bind_rows(arf),root,'arf_channel_contribution',cfg,'n',c('analysis_id','local_cluster'))
 assumptions<-tibble::tibble(item=c('SOFA','urine_output','APACHE','care_preferences','cancer_stage','ICU_free_days','validation_scope'),status=c('modified_creatinine_only_missing_domains_unknown','not_used_check_source_availability','not_validated_do_not_claim_APACHE_comparison','source_orders_do_not_establish_all_treatment_limits','optional_not_inferred_from_advanced_cancer_proxy','not_computed_total_LOS_regression_removed','record_actual_local_or_heldout_design'))
 supplement_write(assumptions,root,'interpretation_limits',cfg)
 invisible(TRUE)
}
run_site_sofa_audit <- function(cfg, root, config_path, baseline=root) {
 output<-file.path(root,'local','supplement_sofa')
 dir.create(file.path(output,'local'),recursive=TRUE,showWarnings=FALSE)
 # Cached runs may still contain a tilde path from an older configuration reader.
 # Use the canonical path for BOTH environment variables and stdout/stderr.
 output<-normalizePath(path.expand(output),winslash='/',mustWork=TRUE)
 # Scientific parameters come from the cached run, not diagnostic defaults.
 envnames<-c('CLIF_SOFA_CONFIG','CLIF_SOFA_BASELINE','CLIF_SOFA_OUTPUT')
 before<-Sys.getenv(envnames,unset=NA_character_)
 on.exit(for(i in seq_along(envnames))if(is.na(before[i]))Sys.unsetenv(envnames[i]) else do.call(Sys.setenv,setNames(list(before[i]),envnames[i])),add=TRUE)
 Sys.setenv(CLIF_SOFA_CONFIG=normalizePath(path.expand(config_path),winslash='/',mustWork=TRUE),CLIF_SOFA_BASELINE=normalizePath(path.expand(baseline),winslash='/',mustWork=TRUE),CLIF_SOFA_OUTPUT=output)
 statuses<-list()
 for(stage in c('audit','extra_source_audit','extend')){
  script<-file.path('diagnostics/sofa_completeness',paste0(stage,'.R'));log<-file.path(output,'local',paste0(stage,'.log'))
  rc<-system2(file.path(R.home('bin'),'Rscript'),c(shQuote(script),shQuote(cfg$site_name)),stdout=log,stderr=log)
  statuses[[stage]]<-tibble::tibble(stage=stage,status=if(rc==0)'completed' else 'failed_inspect_local_log')
  if(rc!=0)break
 }
 status<-dplyr::bind_rows(statuses);supplement_write(status,root,'sofa_source_audit_status',cfg)
 if(!all(status$status=='completed')||nrow(status)!=3)stop('SOFA supplement audit failed; inspect local/supplement_sofa/local logs. No supplement completion manifest was written.')
 x<-read.csv(file.path(output,'reports',paste0(cfg$site_name,'_staged_completeness.csv')))
 # Do not return exact small change counts alongside linked completeness margins.
 x$scores_changed<-x$changed_existing>0;x$changed_existing<-NULL
 x$missing_n<-x$n-x$available
 supplement_write(x,root,'sofa_recovery_components',cfg,c('n','available','missing_n'),'component')
 for(pair in list(c('evidence_availability','hospitalizations'),c('extra_source_counts','unconverted_cv'),c('expanded_lab_window','n'),c('gcs_value_audit','stays'),c('recovery_flags','count'))){
  f<-file.path(output,'reports',paste0(cfg$site_name,'_',pair[1],'.csv'));z<-read.csv(f)
  countcols<-names(z)[vapply(z,is.numeric,logical(1))]
  supplement_write(z,root,paste0('sofa_',pair[1]),cfg,countcols,character())
 }
 # Source vocabulary is returned without patient values or free-text assessments.
 for(kind in c('vasoactive_mapping_flags','vasoactive_units')){
  z<-read.csv(file.path(output,'reports',paste0(cfg$site_name,'_',kind,'.csv')))
  supplement_write(z,root,paste0('sofa_',kind),cfg,intersect(c('rows','stays'),names(z)),character())
 }
 metadata<-read.csv(file.path(output,'reports',paste0(cfg$site_name,'_all_lab_metadata.csv')))
 metadata<-metadata[grepl('po2|fio2|platelet|bilirubin|creatinine',paste(metadata$lab_name,metadata$lab_category),ignore.case=TRUE),,drop=FALSE]
 supplement_write(metadata,root,'sofa_source_lab_catalogue',cfg)
 invisible(TRUE)
}
finalize_site_supplement <- function(root,cfg,phase='derive',baseline=root) {
 src<-file.path(root,'local','supplement_package','share');dest<-file.path(root,'share','supplement')
 if(dir.exists(dest))stop('Supplement already exists; choose a new supplement output directory')
 files<-list.files(src,pattern='csv$',full.names=TRUE)
 stopifnot(length(files)>0)
 for(f in files){z<-readr::read_csv(f,show_col_types=FALSE);stopifnot(!any(names(z)%in%c('patient_id','hospitalization_id','t0','recorded_dttm','admission_dttm','discharge_dttm','event_time','sequence')))}
 dir.create(dest,recursive=TRUE);stopifnot(all(file.copy(files,dest)))
 sources<-c(list.files('diagnostics/sofa_completeness',pattern='[.]R$',full.names=TRUE),list.files('code/lib',pattern='[.]R$',full.names=TRUE))
 jsonlite::write_json(list(schema_version=RUN_SCHEMA,site=cfg$site_name,run_id=cfg$analysis$run_id,phase=phase,status='complete_aggregate_capture_requires_local_review',parameters=cfg$analysis,analysis_source_id=readRDS(file.path(baseline,'local/prepared.rds'))$source_id,supplement_source_md5=as.list(tools::md5sum(sources)),files=as.list(tools::md5sum(list.files(dest,full.names=TRUE))),no_patient_rows=TRUE,source_scores_not_modified=TRUE,privacy='Complementary count masking applied; linked package requires site disclosure review'),file.path(dest,'manifest.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
 invisible(dest)
}
capture_site_supplement <- function(cfg,root=supplement_path(cfg),config_path='config/config.json') {
 baseline<-root
 root<-Sys.getenv('CLIF_SUPPLEMENT_OUTPUT',unset=root)
 if(dir.exists(file.path(root,'share','supplement')))stop('Supplement already exists; do not overwrite a reviewed package')
 prepared<-readRDS(file.path(baseline,'local/prepared.rds'));analyses<-readRDS(Sys.getenv('CLIF_SUPPLEMENT_DERIVATION',unset=file.path(baseline,'local/derivation.rds')))$analyses
 site_supplement_tables(prepared,analyses,cfg,root)
 site_source_and_reporting_tables(prepared,cfg,root)
 run_site_sofa_audit(cfg,root,config_path,baseline)
 ap<-capture_apache_adapted_components(prepared,analyses,cfg,root)
 panel<-descriptive_panel_from_sequence(analyses$h72_primary$fit)
 primary<-analyses$h72_primary
 if(!is.null(primary$consensus_id)){d<-primary$data;d$phenotype<-d$local_cluster;capture_descriptive_tables(prepared,d,cfg,root,panel=panel,apache=ap,phenotypes=paste0('P',seq_len(primary$fit$k)),consensus_id=primary$consensus_id)}else capture_descriptive_tables(prepared,primary$data,cfg,root,panel=panel,apache=ap)
 supplement_write(tibble::tibble(status=if(is.null(primary$consensus_id))'captured_local_baseline_groups_shared_Table2_requires_frozen_model_phase' else 'captured_frozen_shared_phenotypes',path='share/descriptives'),root,'publication_descriptors_status',cfg)
 supplement_readiness(root,cfg)
 finalize_site_supplement(root,cfg,baseline=baseline)
}
site_source_and_reporting_tables <- function(prepared,cfg,root) {
 files<-list.files(cfg$tables_path,recursive=TRUE,full.names=FALSE,pattern='\\.(parquet|csv|fst)$',ignore.case=TRUE)
 bases<-tolower(sub('^clif_','',tools::file_path_sans_ext(basename(files)),ignore.case=TRUE))
 source<-tibble::tibble(table=c('labs','vitals','patient_assessments','respiratory_support','medication_admin_continuous','medication_admin_intermittent','code_status','oncology','patient_procedures','crrt_therapy','intermittent_hemodialysis','intake_output','urine_output'),available=c('labs','vitals','patient_assessments','respiratory_support','medication_admin_continuous','medication_admin_intermittent','code_status','oncology','patient_procedures','crrt_therapy','intermittent_hemodialysis','intake_output','urine_output')%in%bases)
 supplement_write(source,root,'source_table_availability',cfg)
 # Use the site order vocabulary; absence of an order is never absence of preference.
 hit<-files[bases=='code_status'];catalogue<-NULL
 if(length(hit)==1){
  path<-file.path(cfg$tables_path,hit);ext<-tolower(tools::file_ext(path))
  if(ext=='parquet'){
   schema<-clif_parquet_schema(path);cols<-intersect(schema$names,c('code_status_name','code_status_category'))
   if(length(cols))catalogue<-clif_read_parquet(path,columns=cols,distinct=TRUE)
  }else if(!is.null(prepared$tables$code_status))catalogue<-prepared$tables$code_status |> dplyr::select(dplyr::any_of(c('code_status_name','code_status_category'))) |> dplyr::distinct()
 }
 if(!is.null(catalogue)&&ncol(catalogue)){
  catalogue<-with_columns(catalogue,c('code_status_name','code_status_category'))
  cat<-toupper(gsub('[^A-Za-z]','',catalogue$code_status_category));nm<-toupper(catalogue$code_status_name)
  catalogue$interpretation<-ifelse(cat%in%c('DNRDNI','DNARDNI','DNI') & (is.na(nm)|!grepl('DNI|DO NOT INTUBATE',nm)),'canonical_DNI_not_explicit_in_source_name','source_vocabulary_requires_site_confirmation')
  supplement_write(catalogue,root,'code_status_source_catalogue',cfg)
 }else supplement_write(tibble::tibble(status='source_names_not_available_semantics_unverified'),root,'code_status_source_catalogue',cfg)

}
supplement_readiness <- function(root,cfg) {
 source<-file.path(root,'local/supplement_package/share')
 needed<-c('publication_descriptors_status','apache_adapted_components','apache_adapted_provenance','apache_adapted_status','sofa_original_components','sofa_missing_patterns','sofa_recovery_components','sofa_source_audit_status','sofa_evidence_availability','sofa_extra_source_counts','sofa_expanded_lab_window','sofa_gcs_value_audit','sofa_recovery_flags','sofa_vasoactive_mapping_flags','sofa_vasoactive_units','sofa_source_lab_catalogue','oncology_availability','landmark_reconciliation','repeated_admissions','analysis_return_inventory','early_events_by_absorbing_cluster','annual_hourly_evidence','arf_channel_contribution','interpretation_limits','source_table_availability','code_status_source_catalogue')
 inventory<-tibble::tibble(family=needed,phase='local_supplement',status=ifelse(file.exists(file.path(source,paste0(needed,'.csv'))),'captured','missing'))
 # Tiny cohorts may legitimately have no absorbing fit. Report this explicitly.
 inventory$status[inventory$family=='early_events_by_absorbing_cluster' & inventory$status=='missing']<-'not_estimable_check_derivation_status'
 supplement_write(inventory,root,'return_file_inventory',cfg)
 if(any(inventory$status=='missing'))stop('Required site supplement output is missing; no completion manifest written')
 invisible(inventory)
}
# Reproduce diagnostics from local saved data without changing frozen models.
supplement_model_diagnostics <- function(saved,metadata) {
 rows<-list();fits<-list();errors<-list();term_tests<-list()
 for(k in names(saved$cox)){
  fit<-saved$cox[[k]];spec<-if(grepl('^measured_sofa_',k))'measured_sofa' else if(grepl('^unadjusted_',k))'unadjusted' else 'core'
  endpoint<-sub(paste0('^',spec,'_'),'',k)
  meta<-metadata[metadata$model==spec & metadata$endpoint==endpoint,,drop=FALSE]
  n<-if(inherits(fit,'coxph'))fit$n else NA_integer_;events<-if(inherits(fit,'coxph'))fit$nevent else NA_integer_
  p<-NA_real_;delta<-NA_real_;status<-'original_model_not_estimable'
  if(inherits(fit,'coxph')&&nrow(meta)==1&&meta$status=='ok'){
   z<-saved$data;if(spec=='measured_sofa')z<-z[!is.na(z$sofa_total),,drop=FALSE]
   z$event<-cox_event(z,endpoint)
   # Refit only to materialize the original local design; compare all coefficients.
   refit<-tryCatch(suppressWarnings(survival::coxph(formula(fit),data=z,ties='efron',x=TRUE,y=TRUE,model=FALSE)),error=identity)
   if(inherits(refit,'coxph')){
    delta<-max(abs(stats::coef(refit)-stats::coef(fit)))
    same<-identical(names(stats::coef(refit)),names(stats::coef(fit)))&&refit$n==fit$n&&refit$nevent==fit$nevent&&is.finite(delta)&&delta<1e-9
    if(same){
     test<-tryCatch(survival::cox.zph(refit)$table,error=identity)
     if(is.matrix(test)){
      p<-test['GLOBAL','p'];status<-'computed_original_coefficients_reproduced';fits[[k]]<-refit
      term_tests[[k]]<-tibble::tibble(model_key=k,term=rownames(test),chisq=test[,'chisq'],df=test[,'df'],p=test[,'p'],n=n,events_n=events,non_events_n=n-events)
     }
     else {status<-'PH_test_unavailable';errors[[k]]<-conditionMessage(test)}
    }else status<-'diagnostic_reconstruction_mismatch'
   }else {status<-'diagnostic_refit_failed';errors[[k]]<-conditionMessage(refit)}
  }
  rows[[k]]<-tibble::tibble(model_key=k,n=n,events_n=events,non_events_n=n-events,global_PH_p=p,coefficient_max_abs_difference=delta,status=status)
 }
 list(table=dplyr::bind_rows(rows),term_tests=dplyr::bind_rows(term_tests),fits=fits,errors=errors)
}
# Exploratory time-window sensitivity after PH diagnostics; never replaces core HRs.
supplement_piecewise_associations <- function(d,min_cell=10) {
 min_cell<-aggregate_min_cell(min_cell)
 pterms<-grep('^phenotype_P[0-9]+$',names(d),value=TRUE)
 core<-c('age10','male','black','other_or_unknown_race','year','charlson5','advanced','code_nonfull','code_unknown')
 rows<-list();meta<-list();fits<-list()
 for(start in c(0,7))for(endpoint in cox_endpoints()){
  end<-if(start==0)7 else 28;z<-d[d$followup>start,,drop=FALSE]
  z$interval_time<-pmin(z$followup,end)-start
  z$interval_event<-cox_event(z,endpoint,horizon=end)
  terms<-c(pterms,core);terms<-terms[vapply(z[terms],function(x)length(unique(x))>1,logical(1))]
  f<-as.formula(paste('survival::Surv(interval_time,interval_event)~',paste(c(terms,'strata(info_status)','cluster(patient_id)'),collapse='+')))
  warnings<-character();fit<-tryCatch(withCallingHandlers(survival::coxph(f,data=z,ties='efron',x=TRUE,y=TRUE,model=FALSE),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=identity)
  good<-inherits(fit,'coxph')&&!length(warnings)&&all(is.finite(coef(fit)))
  key<-paste(start,end,endpoint,sep='_');fits[[key]]<-fit
  meta[[key]]<-tibble::tibble(interval_start_day=start,interval_end_day=end,endpoint=endpoint,candidate_n=nrow(z),fitted_n=if(inherits(fit,'coxph'))fit$n else 0L,events_n=sum(z$interval_event),non_events_n=sum(!z$interval_event),status=if(good)'ok' else'nonestimable_or_warning',formula=paste(deparse(f),collapse=' '),analysis='exploratory_posthoc_0to7_7to28_days_after_landmark')
  for(term in pterms){
   ph<-sub('^phenotype_','',term)
   safe<-all(vapply(c('P1',ph),function(g){q<-z[z$phenotype==g,,drop=FALSE];sum(q$interval_event)>=min_cell && sum(1-q$interval_event)>=min_cell},logical(1)))
   ok<-good&&safe&&term%in%names(coef(fit))
   rows[[paste(key,term)]]<-tibble::tibble(interval_start_day=start,interval_end_day=end,endpoint=endpoint,phenotype=ph,reference='P1',estimate=if(ok)unname(coef(fit)[term]) else NA_real_,std_error=if(ok)sqrt(vcov(fit)[term,term]) else NA_real_,status=if(!safe)'withheld_small_event_partition' else if(!ok)'nonestimable' else'ok',analysis='exploratory_posthoc_0to7_7to28_days_after_landmark')
  }
 }
 list(associations=dplyr::bind_rows(rows),metadata=dplyr::bind_rows(meta),fits=fits)
}
