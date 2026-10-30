source('code/lib/load.R')
newpath<-Sys.getenv('CLIF_UCMC_REFRESH_PATH',unset='/Users/saborpete/Library/CloudStorage/Box-Box/03-CLIF-2.1/2.1.0')
out<-Sys.getenv('CLIF_UCMC_REFRESH_OUTPUT',unset='output/run_UCMC_etl_refresh_audit_20261026');dir.create(file.path(out,'local'),recursive=TRUE,showWarnings=FALSE);dir.create(file.path(out,'aggregate'),recursive=TRUE,showWarnings=FALSE)
write_json<-function(x,f)jsonlite::write_json(x,file.path(out,f),auto_unbox=TRUE,pretty=TRUE,null='null',digits=16)
write_agg<-function(x,f){stopifnot(!any(names(x)%in%c('patient_id','hospitalization_id','recorded_dttm','t0')));readr::write_csv(x,file.path(out,'aggregate',f))}
p<-readRDS('output/run_UCMC_ccx_pilot_20261021_v3b_joint_primary/local/prepared.rds')
baseline<-readRDS('output/run_UCMC_ccx_pilot_20261021_v3b_joint_primary/local/derivation.rds')$analyses$h72_primary
d<-baseline$data;ids<-p$cohort$hospitalization_id;cfg<-read_run_config();cfg$tables_path<-newpath;cfg$analysis$run_id<-'ucmc_refresh_20261026'
old<-arrow::read_parquet('output/run_UCMC_waterfall_20261022/local/raw_respiratory.parquet');old$recorded_dttm<-utc(old$recorded_dttm)
new<-arrow::open_dataset(file.path(newpath,'clif_respiratory_support.parquet')) |> dplyr::filter(hospitalization_id%in%ids) |> dplyr::collect();new$recorded_dttm<-utc(new$recorded_dttm);new$fio2_set<-standard_fio2(new$fio2_set)
saveRDS(new,file.path(out,'local/new_respiratory_source.rds'))
message('Loaded refreshed source: ',nrow(new),' matched-cohort rows')
# Anchor-table comparison is separate from the respiratory-only paired panel.
tables<-list()
cols<-list(patient=c('patient_id','sex_category','race_category','ethnicity_category','death_dttm'),hospitalization=c('patient_id','hospitalization_id','admission_dttm','discharge_dttm','age_at_admission','discharge_category'),adt=c('hospitalization_id','in_dttm','out_dttm','location_category'),hospital_diagnosis=c('hospitalization_id','diagnosis_code','diagnosis_code_format','poa_present'))
for(nm in names(cols)){
 x<-arrow::open_dataset(file.path(newpath,paste0('clif_',nm,'.parquet'))) |> dplyr::select(dplyr::all_of(cols[[nm]])) |> dplyr::collect()
 for(ts in grep('_dttm$',names(x),value=TRUE))x[[ts]]<-utc(x[[ts]])
 tables[[nm]]<-x
}
built<-build_cohort(tables,cfg);lm<-landmark_cohort(built$cohort,72);eligible<-lm[lm$landmark_reason=='eligible',]
write_json(list(previous_eligible_hospitalizations=nrow(p$cohort),current_eligible_hospitalizations=nrow(built$cohort),same_eligible_ids=setequal(p$cohort$hospitalization_id,built$cohort$hospitalization_id),previous_landmark_hospitalizations=nrow(d),current_landmark_hospitalizations=nrow(eligible),same_landmark_ids=setequal(d$hospitalization_id,eligible$hospitalization_id),new_eligible_n=length(setdiff(built$cohort$hospitalization_id,p$cohort$hospitalization_id)),lost_eligible_n=length(setdiff(p$cohort$hospitalization_id,built$cohort$hospitalization_id))), 'aggregate/cohort_comparison.json')
# Compare source records in the unchanged primary window.
scope<-function(r){r |> dplyr::inner_join(d |> dplyr::select(hospitalization_id,t0),by='hospitalization_id') |> dplyr::mutate(h=floor(hours_between(recorded_dttm,t0)),state=support_state(device_category)) |> dplyr::filter(h>=0,h<72)}
o<-scope(old);n<-scope(new)
for_pair<-function(x,label){x |> dplyr::count(state,name='source_rows') |> dplyr::mutate(version=label)}
write_agg(dplyr::bind_rows(for_pair(o,'previous'),for_pair(n,'refreshed')),'source_state_counts.csv')
# Check whether changed timezone metadata hides a constant timestamp displacement.
old_keys<-unique(paste(o$hospitalization_id,as.numeric(o$recorded_dttm),o$state));timing<-lapply(c(-6,-5,0,5,6),function(shift){key<-unique(paste(n$hospitalization_id,as.numeric(n$recorded_dttm)+shift*3600,n$state));tibble::tibble(new_timestamp_shift_hours=shift,old_keys_matched=sum(old_keys%in%key),old_keys=length(old_keys))})
write_agg(dplyr::bind_rows(timing),'timestamp_alignment.csv')
summary_source<-function(x,label){tibble::tibble(version=label,rows=nrow(x),hospitalizations=dplyr::n_distinct(x$hospitalization_id),missing_device_rows=sum(x$state=='UNKNOWN'),fio2_without_device_rows=sum(x$state=='UNKNOWN'&!is.na(x$fio2_set)),positive_flow_without_device_rows=sum(x$state=='UNKNOWN'&num(x$lpm_set)>0,na.rm=TRUE),known_device_rows=sum(x$state!='UNKNOWN'),room_air_rows=sum(x$state=='ROOM_AIR'),unmapped_populated_device_rows=sum(!is.na(x$device_category)&nzchar(trimws(x$device_category))&x$state=='UNKNOWN'))}
write_agg(dplyr::bind_rows(summary_source(o,'previous'),summary_source(n,'refreshed')),'source_summary.csv')
# All nonrespiratory data and the 650 cohort members stay frozen here.
t<-p$tables;t$respiratory_support<-new
panel<-build_panel(d,t,72,cfg);features<-support_summaries(panel)
sq<-baseline$fit$sequence;oldp<-tibble::tibble(hospitalization_id=rep(baseline$fit$ids,ncol(sq)),h=rep(0:71,each=nrow(sq)),old_joint=as.vector(sq))
oldp$old_support<-sub('\\|.*','',oldp$old_joint);oldp$old_phys<-sub('^[^|]*\\|','',oldp$old_joint)
paired<-panel |> dplyr::left_join(oldp,by=c('hospitalization_id','h'))
stopifnot(nrow(paired)==650*72,!anyNA(paired$old_support),sum(paired$old_support=='UNKNOWN')==14525)
old_unknown<-d$hospitalization_id[d$hours_UNKNOWN==72];new_unknown<-features$hospitalization_id[features$hours_UNKNOWN==72]
write_agg(tibble::tibble(method=c('previous_six_hour_carry','refreshed_six_hour_carry'),n=650,hours=46800,unknown_hours=c(sum(paired$old_support=='UNKNOWN'),sum(paired$support=='UNKNOWN')),unknown_pct=c(100*mean(paired$old_support=='UNKNOWN'),100*mean(paired$support=='UNKNOWN')),entirely_unknown_stays=c(length(old_unknown),length(new_unknown)),arf_positive_stays=c(sum(d$measured_arf),sum(features$measured_arf)),newly_known_old_unknown_hours=c(NA_integer_,sum(paired$old_support=='UNKNOWN'&paired$support!='UNKNOWN')),changed_support_hours=c(NA_integer_,sum(paired$old_support!=paired$support)),changed_physiology_hours=c(NA_integer_,sum(paired$old_phys!=paired$physiology))),'paired_raw_panel_comparison.csv')
write_agg(paired |> dplyr::count(old_support,support,name='hours'),'support_crosswalk.csv')
y<-paired |> dplyr::left_join(d |> dplyr::select(hospitalization_id,admission_year),by='hospitalization_id') |> dplyr::group_by(admission_year) |> dplyr::summarise(hours=dplyr::n(),old_unknown_hours=sum(old_support=='UNKNOWN'),new_unknown_hours=sum(support=='UNKNOWN'),old_unknown_pct=100*mean(old_support=='UNKNOWN'),new_unknown_pct=100*mean(support=='UNKNOWN'),.groups='drop');write_agg(y,'unknown_by_year.csv')
# Presence/absence of direct same-hour source evidence among the old 64 stays.
u<-n[n$hospitalization_id%in%old_unknown,]
write_agg(tibble::tibble(previous_entirely_unknown_n=length(old_unknown),now_any_known_source_device_n=dplyr::n_distinct(u$hospitalization_id[u$state!='UNKNOWN']),now_no_source_rows_n=sum(!old_unknown%in%u$hospitalization_id),now_any_source_fio2_n=dplyr::n_distinct(u$hospitalization_id[!is.na(u$fio2_set)]),remaining_entirely_unknown_n=sum(old_unknown%in%new_unknown)),'previously_entirely_unknown_recovery.csv')
# Prepare one bounded, same-cohort waterfall job for a paired preprocessing check.
bounds<-d |> dplyr::transmute(job='h72_primary',hospitalization_id,trajectory_origin=t0,start=pmax(admission_dttm,t0-6*3600),end=t0+72*3600)
arrow::write_parquet(new,file.path(out,'local/raw_respiratory.parquet'));arrow::write_parquet(bounds,file.path(out,'local/window_bounds.parquet'))
saveRDS(list(data=d,new_raw_panel=panel,paired=paired,features=features,cfg=cfg,old_unknown=old_unknown,baseline_prepared_root='output/run_UCMC_ccx_pilot_20261021_v3b_joint_primary'),file.path(out,'local/paired_audit.rds'))
write_json(list(source_path=newpath,source_respiratory_bytes=file.info(file.path(newpath,'clif_respiratory_support.parquet'))$size,source_respiratory_md5=unname(tools::md5sum(file.path(newpath,'clif_respiratory_support.parquet'))),current_analysis_config_path=read_run_config()$tables_path,new_source_timezone_metadata=arrow::open_dataset(file.path(newpath,'clif_respiratory_support.parquet'))$schema$ToString(),paired_cohort_frozen=TRUE,nonrespiratory_clinical_inputs_frozen=TRUE,existing_results_unchanged=TRUE,new_source_provenance='Values in updated CLIF source table; do not infer bedside observation versus ETL filling without producing-code confirmation',audit_source=as.list(tools::md5sum('diagnostics/etl_refresh/check_ucmc.R'))),'aggregate/audit_manifest.json')
message('Paired refreshed-source audit complete; existing analysis files unchanged')
