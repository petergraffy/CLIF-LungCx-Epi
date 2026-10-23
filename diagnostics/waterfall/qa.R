source("code/lib/load.R")
out<-Sys.getenv("CLIF_WATERFALL_OUTPUT",unset="output/run_UCMC_waterfall_20261022")
result<-readRDS(file.path(out,"local","waterfall_derivation.rds"))
primary<-result$analyses$h72_primary
x<-as.data.frame(arrow::read_parquet(file.path(out,"local","h72_primary_waterfall.parquet")))
raw<-as.data.frame(arrow::read_parquet(file.path(out,"local","raw_respiratory.parquet"))) |>
 dplyr::select(hospitalization_id,recorded_dttm,source_device=device_category,source_flow=lpm_set)
b<-as.data.frame(arrow::read_parquet(file.path(out,"local","window_bounds.parquet"))) |>
 dplyr::filter(job=="h72_primary") |> dplyr::select(hospitalization_id,trajectory_origin)
z<-x |> dplyr::inner_join(b,by="hospitalization_id") |> dplyr::filter(recorded_dttm>=trajectory_origin) |>
 dplyr::left_join(raw,by=c("hospitalization_id","recorded_dttm")) |>
 dplyr::mutate(h=floor(hours_between(recorded_dttm,trajectory_origin)))
conditions<-list(
 positive_flow_with_waterfall_room_air=!is.na(z$source_flow) & z$source_flow>0 & z$device_category=="room air",
 named_source_device_changed=!is.na(z$source_device) & (is.na(z$device_category)|tolower(z$source_device)!=z$device_category),
 filled_fio2_021_on_low_flow=z$device_category=="nasal cannula" & !is.na(z$fio2_set) & z$fio2_set==.21 & !z$raw_fio2_observed_at_timestamp,
 raw_positive_flow_with_fio2_021=!is.na(z$source_flow) & z$source_flow>0 & !is.na(z$fio2_set) & z$fio2_set==.21)
checks<-dplyr::bind_rows(lapply(names(conditions),function(label) {
 selected_rows<-z[which(conditions[[label]]),]
 saveRDS(selected_rows,file.path(out,"local",paste0(label,"_review.rds")))
 tibble::tibble(check=label,rows=nrow(selected_rows),hours=nrow(dplyr::distinct(selected_rows,hospitalization_id,h)),hospitalizations=dplyr::n_distinct(selected_rows$hospitalization_id))
}))
readr::write_csv(checks,file.path(out,"aggregate","waterfall_source_conflict_checks.csv"));print(checks)
transitions<-primary$panel |> dplyr::filter(!terminal,physiology!=old_physiology) |>
 dplyr::count(old_physiology,physiology,name="hours",sort=TRUE)
readr::write_csv(transitions,file.path(out,"aggregate","h72_primary_physiology_changes.csv"));print(transitions)
for(id in names(result$analyses)) {
 a<-result$analyses[[id]]
 stopifnot(nrow(a$panel)==nrow(a$data)*a$window,setequal(a$fit$ids,a$data$hospitalization_id),!anyNA(a$fit$labels))
 if(a$variant!="absorbing")stopifnot(all(a$data$time_from_landmark>0),!any(a$panel$terminal))
}
input<-jsonlite::fromJSON(file.path(out,"local","input_manifest.json"))
baseline<-readRDS(file.path(input$baseline_root,"local","derivation.rds"))$analyses$h72_primary
stopifnot(setequal(primary$data$hospitalization_id,baseline$data$hospitalization_id),
 identical(sort(table(primary$data$event_type)),sort(table(baseline$data$event_type))))
old_unknown <- baseline$data$hospitalization_id[baseline$data$hours_UNKNOWN==72]
new_unknown <- primary$data$hospitalization_id[primary$data$hours_UNKNOWN==72]
jsonlite::write_json(list(status="passed",all_analyses_complete=length(result$analyses),
 primary_n=nrow(primary$data),primary_events_unchanged=TRUE,
 entirely_unknown_encounters_unchanged=setequal(old_unknown,new_unknown),
 entirely_unknown_before=length(old_unknown),entirely_unknown_after=length(new_unknown),numeric_bfill=FALSE,
 source_md5=as.list(tools::md5sum(list.files("diagnostics/waterfall",full.names=TRUE,pattern="\\.(R|py)$"))),
 synthetic_boundary_tests="test_boundaries.py",prognostic_validation="not_repeated_with_waterfall_inside_training_folds"),
 file.path(out,"local","comparison_qa.json"),auto_unbox=TRUE,pretty=TRUE)
message("Real-data cohort, event, sequence and boundary QA passed")
