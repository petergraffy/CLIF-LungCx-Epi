# Run from repository root. Reads the configured completed run; does not alter it.
# All source rows, identifiers, and exploratory assignments stay in ignored output.
source("code/lib/load.R")
cfg <- read_run_config()
root <- init_run(cfg)
out <- Sys.getenv("CLIF_SUPPORT_AUDIT_OUTPUT", unset = "output/run_UCMC_support_audit_20261022")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
run <- readRDS(file.path(root,"local","derivation.rds"))
if (!identical(run$source_id, code_fingerprint())) stop("Audit code differs from baseline derivation; use its matching source version")
ids <- readRDS(file.path(root,"local","prepared.rds"))$cohort$hospitalization_id
files <- list.files(cfg$tables_path,recursive=TRUE,full.names=TRUE,pattern="\\.(parquet|csv|fst)$")
base <- tolower(sub("^clif_","",tools::file_path_sans_ext(basename(files)),ignore.case=TRUE))
path <- files[base=="respiratory_support"]
if(length(path)!=1L)stop("Expected exactly one source respiratory_support table")
ext <- tolower(tools::file_ext(path))
if(ext=="parquet") {
  source_rows <- arrow::open_dataset(path) |>
    dplyr::filter(hospitalization_id %in% ids) |> dplyr::collect()
} else {
  source_rows <- if(ext=="fst")fst::read_fst(path) else readr::read_csv(path,show_col_types=FALSE)
  source_rows <- source_rows[source_rows$hospitalization_id %in% ids,]
}
source_rows$recorded_dttm <- utc(source_rows$recorded_dttm)
source_rows <- with_columns(source_rows,c("device_name","mode_category","artificial_airway","peep_set","fio2_set","lpm_set"))
saveRDS(source_rows,file.path(out,"respiratory_source_local.rds"))
for(script in c("diagnose.R","oxygen_evidence.R","exploratory_oxygen_state.R","plot.R"))
  source(file.path("diagnostics","support_missingness",script))
baseline <- run$analyses$h72_primary$data
hours <- readRDS(file.path(out,"hour_audit_local.rds"))
checks <- readr::read_csv(file.path(out,"oxygen_evidence_summary.csv"),show_col_types=FALSE)
stopifnot(nrow(hours)==nrow(baseline)*72L,
 sum(hours$support=="UNKNOWN")==sum(baseline$hours_UNKNOWN),
 sum(checks$observed_oxygen_hours+checks$carried_oxygen_hours+checks$residual_unknown_hours)==sum(checks$unknown_hours))
jsonlite::write_json(list(baseline_run=cfg$analysis$run_id,site=cfg$site_name,
 baseline_source_id=run$source_id,source_file=path,source_bytes=file.info(path)$size,
 audit_code_md5=as.list(tools::md5sum(list.files("diagnostics",recursive=TRUE,full.names=TRUE,pattern="\\.R$"))),
 exploratory_only=TRUE,primary_unchanged=TRUE),file.path(out,"audit_provenance.json"),pretty=TRUE,auto_unbox=TRUE)
message("Support audit complete. All outputs remain local: ",out)
