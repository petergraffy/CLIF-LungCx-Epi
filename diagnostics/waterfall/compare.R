source("code/lib/load.R")
cfg<-read_run_config(Sys.getenv("CLIF_RUN_CONFIG",unset="config/config.json"));out<-Sys.getenv("CLIF_WATERFALL_OUTPUT",unset="output/run_UCMC_waterfall_20261022")
input<-jsonlite::fromJSON(file.path(out,"local","input_manifest.json"))
baseline<-readRDS(file.path(input$baseline_root,"local","derivation.rds"))
prepared<-readRDS(file.path(input$baseline_root,"local","prepared.rds"))
stopifnot(identical(baseline$source_id,code_fingerprint()))
for(folder in c("local","aggregate","figures"))dir.create(file.path(out,folder),recursive=TRUE,showWarnings=FALSE)
parts<-function(x,i) vapply(strsplit(x,"|",fixed=TRUE),function(z)if(length(z)>=i)z[i] else "HU_CU",character(1))
old_panel<-function(a) {
 sq<-a$fit$sequence
 tibble::tibble(hospitalization_id=rep(a$fit$ids,ncol(sq)),h=rep(seq_len(ncol(sq))-1L,each=nrow(sq)),
  old_joint=as.vector(sq),old_support=parts(as.vector(sq),1),old_physiology=parts(as.vector(sq),2))
}
plans<-tibble::tibble(analysis_id=c(names(baseline$analyses),"h72_waterfall_support_with_raw_arf"),
 baseline_id=c(names(baseline$analyses),"h72_primary"))
results<-list();summaries<-list();year_summaries<-list();crosswalks<-list();clinical_summaries<-list();stability<-NULL
for(i in seq_len(nrow(plans))) {
 id<-plans$analysis_id[i];base_id<-plans$baseline_id[i];a<-baseline$analyses[[base_id]]
 message("Deriving waterfall ",id)
 job<-if(base_id %in% c("h72_absorbing","h24_landmark_sensitivity","h48_landmark_sensitivity","h72_arf_onset"))base_id else "h72_primary"
 rs<-as.data.frame(arrow::read_parquet(file.path(out,"local",paste0(job,"_waterfall.parquet"))))
 rs$recorded_dttm<-utc(rs$recorded_dttm)
 rs<-rs[rs$hospitalization_id %in% a$data$hospitalization_id,]
 limits<-a$data |> dplyr::transmute(hospitalization_id,
  start=pmax(admission_dttm,t0-cfg$analysis$support_carry_hours*3600),end=pmin(t0+a$window*3600,event_time))
 bounded<-rs |> dplyr::left_join(limits,by="hospitalization_id")
 stopifnot(all(bounded$recorded_dttm>=bounded$start & bounded$recorded_dttm<bounded$end))
 tables<-prepared$tables;tables$respiratory_support<-rs
 panel<-build_panel(a$data,tables,a$window,cfg,absorbing=a$variant=="absorbing",
  estimate_low_flow=a$variant=="low_flow_fio2",carry=switch(a$variant,arf_carry_0=0,arf_carry_12=12,NULL))
 comparison_a<-if(a$fit$field=="support")baseline$analyses$h72_primary else a
 panel<-panel |> dplyr::left_join(old_panel(comparison_a),by=c("hospitalization_id","h"))
 if(id=="h72_waterfall_support_with_raw_arf") {
  panel$physiology<-panel$old_physiology
  panel$hypox<-vapply(strsplit(panel$physiology,"_",fixed=TRUE),function(z)z[1],character(1))
  panel$hyper<-vapply(strsplit(panel$physiology,"_",fixed=TRUE),function(z)z[2],character(1))
  panel$joint<-paste(panel$support,panel$physiology,sep="|")
 }
 # Processed/scaffold rows are not raw observations. Track actual matching
 # device documentation separately from all inferred or filled hourly states.
 direct<-rs |> dplyr::filter(raw_device_observed_at_timestamp) |>
  dplyr::left_join(a$data |> dplyr::select(hospitalization_id,t0),by="hospitalization_id") |>
  dplyr::mutate(h=floor(hours_between(recorded_dttm,t0)),support=support_state(device_category))
 direct_keys<-paste(direct$hospitalization_id,direct$h,direct$support)
 panel$support_observed<-paste(panel$hospitalization_id,panel$h,panel$support) %in% direct_keys & !panel$terminal
 panel$support_carried<-panel$support!="UNKNOWN" & !panel$support_observed & !panel$terminal
 stopifnot(nrow(panel)==nrow(a$data)*a$window,!anyNA(panel$old_support))
 if(a$variant!="absorbing")stopifnot(!any(panel$support %in% TERMINAL_STATES),all(a$data$time_from_landmark>0))
 if(a$variant=="absorbing") {
  expected<-a$data |> dplyr::select(hospitalization_id,t0,event_time,event_type)
  check<-panel |> dplyr::left_join(expected,by="hospitalization_id")
  terminal<-floor(hours_between(check$event_time,check$t0))<=check$h
  correct<-c(death="DEATH",hospice="HOSPICE",live_discharge="LIVE_DISCHARGE")[check$event_type]
  stopifnot(all(check$support[terminal]==correct[terminal]))
 }
 features<-support_summaries(panel)
 d<-a$data |> dplyr::select(-dplyr::all_of(c("local_cluster",intersect(setdiff(names(features),"hospitalization_id"),names(a$data))))) |>
  dplyr::left_join(features,by="hospitalization_id")
 field<-if(a$variant=="support_only")"support" else "joint"
 fit<-sequence_fit(panel,field,cfg$analysis$k_max)
 stopifnot(setequal(fit$ids,a$fit$ids),!anyDuplicated(fit$ids))
 d$local_cluster<-fit$labels[match(d$hospitalization_id,fit$ids)]
 assignments<-d |> dplyr::select(hospitalization_id,local_cluster)
 crosswalks[[id]]<-a$data |> dplyr::select(hospitalization_id,original_cluster=local_cluster) |>
  dplyr::left_join(assignments,by="hospitalization_id") |> dplyr::count(original_cluster,local_cluster,name="n") |>
  dplyr::mutate(analysis_id=id)
 compared_labels<-fit$labels[match(a$fit$ids,fit$ids)]
 summaries[[id]]<-tibble::tibble(analysis_id=id,baseline_analysis=base_id,n=nrow(d),
  original_k=a$fit$k,waterfall_k=fit$k,original_unknown_hours=sum(panel$old_support=="UNKNOWN"),
  waterfall_unknown_hours=sum(panel$support=="UNKNOWN"),total_hours=nrow(panel),
  newly_known_hours=sum(panel$old_support=="UNKNOWN" & panel$support!="UNKNOWN"),
  newly_unknown_hours=sum(panel$old_support!="UNKNOWN" & panel$support=="UNKNOWN"),
  raw_device_observed_hours=sum(panel$support_observed),filled_or_inferred_support_hours=sum(panel$support_carried),
  physiology_changed_hours=sum(!panel$terminal & panel$physiology!=panel$old_physiology),
  positive_arf_hospitalizations=sum(d$measured_arf),
  original_silhouette=max(a$fit$silhouette$avg_silhouette),waterfall_silhouette=max(fit$silhouette$avg_silhouette),
  adjusted_rand=adjusted_rand(a$fit$labels,compared_labels),
  post_landmark_deaths=if(a$variant=="absorbing")NA_integer_ else sum(d$event_type=="death"))
 year_summaries[[id]]<-panel |> dplyr::left_join(d |> dplyr::select(hospitalization_id,admission_year),by="hospitalization_id") |>
  dplyr::group_by(admission_year) |> dplyr::summarise(hours=dplyr::n(),
   original_unknown_hours=sum(old_support=="UNKNOWN"),waterfall_unknown_hours=sum(support=="UNKNOWN"),.groups="drop") |>
  dplyr::mutate(analysis_id=id)
 clinical_summaries[[id]]<-d |> dplyr::group_by(local_cluster) |> dplyr::summarise(n=dplyr::n(),
  hospital_death_n=sum(event_type=="death"),hospice_n=sum(event_type=="hospice"),
  other_live_discharge_n=sum(event_type=="live_discharge"),arf_positive_n=sum(measured_arf),
  unknown_support_pct=100*mean(hours_UNKNOWN/a$window),measured_sofa_n=sum(!is.na(sofa_total)),.groups="drop") |>
  dplyr::mutate(analysis_id=id)
 readr::write_csv(fit$silhouette,file.path(out,"aggregate",paste0(id,"_silhouette.csv")))
 if(id %in% c("h72_primary","h72_waterfall_support_with_raw_arf")) {
  stability<-bootstrap_stability(fit,d$patient_id[match(fit$ids,d$hospitalization_id)],cfg$analysis$bootstrap_reps,cfg$analysis$seed)
  stability_summary<-stability |> dplyr::group_by(local_cluster) |> dplyr::summarise(
   completed_replicates=dplyr::n(),mean_jaccard=mean(jaccard),lower_jaccard=quantile(jaccard,.025),
   mean_adjusted_rand=mean(adjusted_rand),.groups="drop")
  readr::write_csv(stability_summary,file.path(out,"aggregate",paste0(id,"_stability.csv")))
 }
 if(id=="h72_primary") {
  burden<-panel |> dplyr::inner_join(assignments,by="hospitalization_id") |>
   dplyr::mutate(rank=match(support,SUPPORT_STATES)-1L,rank=ifelse(support %in% c("OTHER","UNKNOWN"),NA_real_,rank)) |>
   dplyr::group_by(local_cluster) |> dplyr::summarise(burden=mean(rank,na.rm=TRUE),.groups="drop") |>
   dplyr::arrange(burden,local_cluster)
  d$phenotype<-paste0("C",d$local_cluster);reference<-paste0("C",burden$local_cluster[1])
  for(severity in c("measured_components","missing_zero_sensitivity")) {
   x<-d;if(severity=="missing_zero_sensitivity")x$sofa_total<-x$sofa_missing_zero_sensitivity
   model<-fit_associations(x,reference)
   for(kind in names(model))readr::write_csv(model[[kind]],file.path(out,"aggregate",paste0(id,"_",severity,"_",kind,".csv")))
  }
  readr::write_csv(aj_grid(d),file.path(out,"aggregate","h72_primary_cumulative_incidence.csv"))
 }
 fit$distance<-NULL
 results[[id]]<-list(data=d,fit=fit,panel=panel,window=a$window,variant=a$variant,
  processing="clifpy_waterfall_bfill_false_bounded",baseline_analysis=base_id)
 saveRDS(list(analyses=results,comparison=dplyr::bind_rows(summaries)),file.path(out,"local","comparison_progress.rds"))
}
summary<-dplyr::bind_rows(summaries) |> dplyr::mutate(original_unknown_pct=100*original_unknown_hours/total_hours,
 waterfall_unknown_pct=100*waterfall_unknown_hours/total_hours)
readr::write_csv(summary,file.path(out,"aggregate","analysis_comparison.csv"))
readr::write_csv(dplyr::bind_rows(year_summaries),file.path(out,"aggregate","unknown_by_year.csv"))
readr::write_csv(dplyr::bind_rows(crosswalks),file.path(out,"aggregate","cluster_crosswalk.csv"))
readr::write_csv(dplyr::bind_rows(clinical_summaries),file.path(out,"aggregate","cluster_outcomes.csv"))
saveRDS(list(analyses=results,comparison=summary,baseline_source_id=baseline$source_id),file.path(out,"local","waterfall_derivation.rds"))
print(summary |> dplyr::select(analysis_id,n,original_k,waterfall_k,original_unknown_pct,waterfall_unknown_pct,adjusted_rand),n=20)
message("Waterfall comparison complete; original baseline remains unchanged")
