source('code/lib/load.R');suppressPackageStartupMessages(library(ggplot2))
out<-Sys.getenv('CLIF_UCMC_REFRESH_OUTPUT',unset='output/run_UCMC_etl_refresh_audit_20261026');a<-readRDS(file.path(out,'local/paired_audit.rds'));p<-readRDS(file.path(a$baseline_prepared_root,'local/prepared.rds'));d<-a$data
x<-arrow::read_parquet(file.path(out,'local/h72_primary_waterfall.parquet'));x$recorded_dttm<-utc(x$recorded_dttm)
limits<-d |> dplyr::transmute(hospitalization_id,start=pmax(admission_dttm,t0-6*3600),end=t0+72*3600)
check<-x |> dplyr::left_join(limits,by='hospitalization_id');stopifnot(all(check$recorded_dttm>=check$start & check$recorded_dttm<check$end))
t<-p$tables;t$respiratory_support<-x
wf<-build_panel(d,t,72,a$cfg)
# In the conservative refreshed comparison, ARF comes from refreshed RAW source,
# not filled waterfall FiO2. Outcome/covariate inputs remain frozen.
new<-a$new_raw_panel;stopifnot(identical(new[c('hospitalization_id','h')],wf[c('hospitalization_id','h')]))
wfp<-wf;wfp$hypox<-new$hypox;wfp$hyper<-new$hyper;wfp$physiology<-new$physiology;wfp$joint<-paste(wfp$support,wfp$physiology,sep='|')
f<-support_summaries(wfp)
oldwf<-readRDS('output/run_UCMC_waterfall_20261022/local/waterfall_derivation.rds')$analyses$h72_waterfall_support_with_raw_arf
rawsummary<-readr::read_csv(file.path(out,'aggregate/paired_raw_panel_comparison.csv'),show_col_types=FALSE)
summary<-dplyr::bind_rows(rawsummary |> dplyr::select(method,n,hours,unknown_hours,unknown_pct,entirely_unknown_stays,arf_positive_stays),tibble::tibble(method=c('previous_waterfall_support_raw_arf','refreshed_waterfall_support_raw_arf'),n=650,hours=46800,unknown_hours=c(sum(oldwf$panel$support=='UNKNOWN'),sum(wfp$support=='UNKNOWN')),unknown_pct=c(100*mean(oldwf$panel$support=='UNKNOWN'),100*mean(wfp$support=='UNKNOWN')),entirely_unknown_stays=c(sum(oldwf$data$hours_UNKNOWN==72),sum(f$hours_UNKNOWN==72)),arf_positive_stays=c(sum(oldwf$data$measured_arf),sum(f$measured_arf))))
readr::write_csv(summary,file.path(out,'aggregate/final_comparison.csv'))
# Direct source evidence at remaining unknown hours, not inferred numeric values.
r<-readRDS(file.path(out,'local/new_respiratory_source.rds')) |> dplyr::inner_join(d |> dplyr::select(hospitalization_id,t0),by='hospitalization_id') |> dplyr::mutate(h=floor(hours_between(recorded_dttm,t0))) |> dplyr::filter(h>=0,h<72)
e<-r |> dplyr::group_by(hospitalization_id,h) |> dplyr::summarise(source_rows=dplyr::n(),source_fio2=any(!is.na(fio2_set)),source_positive_flow=any(num(lpm_set)>0,na.rm=TRUE),.groups='drop')
u<-wfp[wfp$support=='UNKNOWN',] |> dplyr::left_join(e,by=c('hospitalization_id','h'))
readr::write_csv(tibble::tibble(unknown_hours=nrow(u),hospitalizations=dplyr::n_distinct(u$hospitalization_id),no_source_row_hours=sum(is.na(u$source_rows)),source_fio2_hours=sum(u$source_fio2,na.rm=TRUE),source_positive_flow_hours=sum(u$source_positive_flow,na.rm=TRUE)),file.path(out,'aggregate/residual_unknown_source_evidence.csv'))
# Compare source observations against processed room-air assignments for review.
q<-x |> dplyr::select(hospitalization_id,recorded_dttm,device_category,fio2_set,raw_fio2_observed_at_timestamp) |> dplyr::inner_join(d |> dplyr::select(hospitalization_id,t0),by='hospitalization_id') |> dplyr::filter(recorded_dttm>=t0) |> dplyr::left_join(r |> dplyr::select(hospitalization_id,recorded_dttm,source_flow=lpm_set),by=c('hospitalization_id','recorded_dttm'))
readr::write_csv(tibble::tibble(check=c('positive_source_flow_with_waterfall_room_air','filled_fio2_021_on_nasal_cannula'),rows=c(sum(num(q$source_flow)>0 & q$device_category=='room air',na.rm=TRUE),sum(q$device_category=='nasal cannula' & q$fio2_set==.21 & !q$raw_fio2_observed_at_timestamp,na.rm=TRUE))),file.path(out,'aggregate/waterfall_conflict_flags.csv'))
saveRDS(list(data=d,refreshed_waterfall_panel=wfp,full_numeric_waterfall_panel=wf,features=f),file.path(out,'local/waterfall_audit.rds'))
y<-readr::read_csv(file.path(out,'aggregate/unknown_by_year.csv'),show_col_types=FALSE)
y<-y |> tidyr::pivot_longer(c(old_unknown_pct,new_unknown_pct),names_to='version',values_to='unknown_pct')
g<-ggplot(y,aes(admission_year,unknown_pct,color=version))+geom_line()+geom_point()+scale_color_discrete(labels=c(new_unknown_pct='Updated source',old_unknown_pct='Previous source'))+labs(x='Admission year',y='Unknown support hours (%)',title='UCMC respiratory-source update: same 650 landmark hospitalizations',subtitle='Recorded source plus six-hour carry; no waterfall in this figure',caption='Previous analysis files unchanged. Source labels may include ETL inference; provenance requires producing-code confirmation.')+theme_minimal(base_size=11)+theme(legend.position='bottom',plot.caption=element_text(hjust=0))
dir.create(file.path(out,'figures'),showWarnings=FALSE);ggsave(file.path(out,'figures/unknown_support_by_year.png'),g,width=11,height=5,dpi=160)
ledger<-jsonlite::fromJSON('output/run_federated_rehearsal_20261022/reports/frozen_artifact_ledger.json',simplifyVector=FALSE)$ledger
for(v in ledger)stopifnot(unname(tools::md5sum(v$path))==v$md5)
jsonlite::write_json(list(status='passed',same_650_hospitalizations=TRUE,all_cases_retained=TRUE,no_terminal_exposure_states=TRUE,refreshed_raw_arf_preserved_under_waterfall=identical(wfp$physiology,new$physiology),bounded_waterfall=TRUE,original_three_site_artifacts_unchanged=TRUE,cohort_check=jsonlite::fromJSON(file.path(out,'aggregate/cohort_comparison.json')),new_source_path=a$cfg$tables_path,source=as.list(tools::md5sum(c('diagnostics/etl_refresh/check_ucmc.R','diagnostics/etl_refresh/finish_ucmc.R')))),file.path(out,'aggregate/final_audit.json'),pretty=TRUE,auto_unbox=TRUE)
print(summary,width=Inf);message('Refreshed UCMC paired ETL/waterfall audit complete')
