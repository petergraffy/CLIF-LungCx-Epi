source("code/lib/load.R")
cfg <- read_run_config(Sys.getenv("CLIF_RUN_CONFIG",unset="config/config.json"))
baseline_root <- init_run(cfg)
out <- Sys.getenv("CLIF_WATERFALL_OUTPUT",unset="output/run_UCMC_waterfall_20261022")
dir.create(file.path(out,"local"),recursive=TRUE,showWarnings=FALSE)
baseline <- readRDS(file.path(baseline_root,"local","derivation.rds"))
prepared <- readRDS(file.path(baseline_root,"local","prepared.rds"))
stopifnot(identical(baseline$source_id,code_fingerprint()))
jobs <- c("h72_primary","h72_absorbing","h24_landmark_sensitivity","h48_landmark_sensitivity","h72_arf_onset")
bounds <- dplyr::bind_rows(lapply(jobs,function(id) {
 a<-baseline$analyses[[id]];d<-a$data
 tibble::tibble(job=id,hospitalization_id=d$hospitalization_id,
  trajectory_origin=d$t0,
  start=pmax(d$admission_dttm,d$t0-cfg$analysis$support_carry_hours*3600),
  end=pmin(d$t0+a$window*3600,d$event_time))
}))
stopifnot(all(bounds$end>bounds$start),!anyDuplicated(bounds[c("job","hospitalization_id")]))
files<-list.files(cfg$tables_path,recursive=TRUE,full.names=TRUE,pattern="\\.parquet$")
base<-tolower(sub("^clif_","",tools::file_path_sans_ext(basename(files)),ignore.case=TRUE))
path<-files[base=="respiratory_support"];stopifnot(length(path)==1L)
ids<-prepared$cohort$hospitalization_id
dataset<-arrow::open_dataset(path)
raw<-dataset |> dplyr::filter(hospitalization_id %in% ids) |> dplyr::collect()
raw$recorded_dttm<-source_time(raw$recorded_dttm,dataset$schema$GetFieldByName("recorded_dttm")$type)
# Retain the primary analysis's row-wise FiO2 normalization and bounds.
raw$fio2_set<-standard_fio2(raw$fio2_set)
arrow::write_parquet(raw,file.path(out,"local","raw_respiratory.parquet"))
arrow::write_parquet(bounds,file.path(out,"local","window_bounds.parquet"))
jsonlite::write_json(list(baseline_run=cfg$analysis$run_id,baseline_root=baseline_root,
 baseline_source_id=baseline$source_id,site=cfg$site_name,jobs=jobs,
 pre_window_lookback_hours=cfg$analysis$support_carry_hours,
 boundary_rule="same-admission six-hour lookback; no input or output at/after landmark or terminal event",
 numeric_bfill=FALSE,categorical_within_window_lookahead=TRUE,
 original_clinical_covariates_frozen=TRUE,source_file=path),
 file.path(out,"local","input_manifest.json"),auto_unbox=TRUE,pretty=TRUE)
message("Prepared five bounded waterfall jobs without changing the baseline run")
