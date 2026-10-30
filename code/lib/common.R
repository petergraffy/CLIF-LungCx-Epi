# Shared contracts for the new CCX run. Never source archived analysis code.
RUN_SCHEMA <- "ccx-2026-10-21-v3"
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

run_defaults <- function() list(
  start_date = "2018-01-01", end_date = "2025-01-01", seed = 20261021L,
  windows = c(24L, 48L, 72L), support_carry_hours = 6,
  physiology_carry_hours = 6, pairing_hours = 1, hyper_pairing_hours = 2,
  bootstrap_reps = 100L,
  validation_folds = 5L, validation_bootstrap_reps = 200L,
  min_cell = 0L, k_max = 8L, prediction_days = 28,
  run_id = "ccx_revision_v3", consensus_k = NULL,
  code_status_lookback_hours = 0, sf_threshold = 315,
  consensus_min_overlap = 0.5
)

analysis_specifications <- function(cfg) {
  early <- sort(setdiff(cfg$analysis$windows, 72L))
  tibble::tibble(
    analysis_id = c("h72_primary", "h72_absorbing", "h72_support_only", "h72_arf_onset",
      paste0("h", early, "_landmark_sensitivity"), "h72_low_flow_fio2", "h72_arf_carry_0", "h72_arf_carry_12"),
    window = c(rep(72L, 4), early, rep(72L, 3)),
    variant = c("primary", "absorbing", "support_only", "arf_onset",
      rep("earlier_landmark", length(early)), "low_flow_fio2", "arf_carry_0", "arf_carry_12"),
    analysis_role = c("primary", rep("sensitivity", 6 + length(early))),
    trajectory_origin = c(rep("icu_admission", 3), "first_detected_arf", rep("icu_admission", 3 + length(early))),
    clustering_inputs = c("joint_support_arf", "joint_support_arf_terminal_events", "support_only",
      rep("joint_support_arf", 4 + length(early))),
    outcome_modeling = c("post_landmark", "descriptive_only", rep("post_landmark", 5 + length(early))))
}

read_run_config <- function(path = "config/config.json") {
  if (!file.exists(path)) stop("Create config/config.json from config/config_template.json")
  cfg <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  for (nm in c("site_name", "tables_path", "repo"))
    if (is.null(cfg[[nm]]) || !nzchar(cfg[[nm]])) stop("Missing config field: ", nm)
  cfg$analysis <- modifyList(run_defaults(), cfg$analysis %||% list())
  cfg$analysis$min_cell <- aggregate_min_cell(cfg$analysis$min_cell)
  # Retired pilot parameter: support coverage is descriptive, never eligibility.
  cfg$analysis$minimum_support_coverage <- NULL
  a <- cfg$analysis
  study_bounds <- utc(c(a$start_date, a$end_date))
  if (length(study_bounds) != 2L || anyNA(study_bounds) || study_bounds[1] >= study_bounds[2])
    stop("Study start/end must be valid dates or timestamps, with start before end")
  stopifnot(all(a$windows %in% c(24L, 48L, 72L)), 72L %in% a$windows,
            a$support_carry_hours >= 0, a$physiology_carry_hours >= 0,
            a$min_cell >= 0, a$bootstrap_reps >= 0, a$validation_folds >= 2,
            a$validation_bootstrap_reps >= 2, a$k_max >= 2,
            a$prediction_days > 0)
  cfg$site_name <- gsub("[^A-Za-z0-9_-]", "_", cfg$site_name)
  cfg$analysis$run_id <- gsub("[^A-Za-z0-9_-]", "_", a$run_id)
  cfg
}

utc <- function(x) {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = "UTC"))
  if (is.numeric(x)) return(as.POSIXct(x, origin = "1970-01-01", tz = "UTC"))
  # ymd_hms(truncated=2) still requires an hour. Study boundaries are
  # commonly date-only; parse them explicitly rather than silently yielding NA.
  if (is.character(x)) {
    date_only <- !is.na(x) & grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", trimws(x))
    x[date_only] <- paste0(trimws(x[date_only]), " 00:00:00")
  }
  suppressWarnings(lubridate::ymd_hms(x, quiet = TRUE, truncated = 2, tz = "UTC"))
}
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
hours_between <- function(x, y) as.numeric(difftime(x, y, units = "hours"))
need_columns <- function(x, required, table) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(table, " missing required columns: ", paste(missing, collapse = ", "))
  x
}
with_columns <- function(x, cols) {
  for (nm in setdiff(cols, names(x))) x[[nm]] <- NA
  x
}

load_clif <- function(cfg) {
  required <- list(
    patient = c("patient_id", "sex_category", "race_category"),
    hospitalization = c("patient_id", "hospitalization_id", "admission_dttm", "discharge_dttm", "discharge_category", "age_at_admission"),
    adt = c("hospitalization_id", "in_dttm", "out_dttm", "location_category"),
    hospital_diagnosis = c("hospitalization_id", "diagnosis_code", "diagnosis_code_format", "poa_present"),
    respiratory_support = c("hospitalization_id", "recorded_dttm", "device_category"),
    vitals = c("hospitalization_id", "recorded_dttm", "vital_category", "vital_value"),
    labs = c("hospitalization_id", "lab_result_dttm", "lab_category", "lab_value_numeric"),
    medication_admin_continuous = c("hospitalization_id", "admin_dttm", "med_category", "med_dose", "med_dose_unit"),
    patient_assessments = c("hospitalization_id", "recorded_dttm", "assessment_category", "numerical_value")
  )
  files <- list.files(cfg$tables_path, recursive = TRUE, full.names = TRUE,
                      pattern = "\\.(csv|parquet|fst)$", ignore.case = TRUE)
  bases <- tolower(sub("^clif_", "", tools::file_path_sans_ext(basename(files)), ignore.case = TRUE))
  out <- list()
  cohort_ids <- patient_ids <- NULL
  optional_columns <- list(patient = c("death_dttm", "ethnicity_category"),
    respiratory_support = c("fio2_set", "lpm_set", "flow_rate"),
    code_status = c("patient_id", "start_dttm", "code_status_category"),
    oncology = c("hospitalization_id", "histology", "stage", "ecog", "therapy_line", "therapy_intent", "molecular_subtype", "treatment_type"))
  for (nm in c(names(required), "code_status", "oncology")) {
    if (nm == "respiratory_support") {
      eligible <- build_cohort(out, cfg)$cohort
      cohort_ids <- eligible$hospitalization_id; patient_ids <- unique(eligible$patient_id)
    }
    hits <- files[bases == nm]
    if (!length(hits) && !nm %in% names(required)) next
    if (length(hits) != 1) stop("Expected one CLIF table for ", nm, "; found ", length(hits))
    ext <- tolower(tools::file_ext(hits))
    if (ext == "parquet") {
      query <- arrow::open_dataset(hits)
      columns <- query$schema$names
      selected <- columns[tolower(columns) %in% c(required[[nm]], optional_columns[[nm]])]
      query <- dplyr::select(query, dplyr::all_of(selected))
      # Push cohort filtering into Arrow before materializing large clinical tables.
      hid <- selected[tolower(selected) == "hospitalization_id"]
      pid <- selected[tolower(selected) == "patient_id"]
      if (!is.null(cohort_ids) && length(hid)) query <- dplyr::filter(query, .data[[hid]] %in% cohort_ids)
      else if (!is.null(patient_ids) && length(pid)) query <- dplyr::filter(query, .data[[pid]] %in% patient_ids)
      x <- as.data.frame(dplyr::collect(query))
    } else x <- switch(ext,
      csv = readr::read_csv(hits, show_col_types = FALSE, col_types = readr::cols(.default = readr::col_character())),
      fst = fst::read_fst(hits))
    names(x) <- tolower(names(x))
    if (nm %in% names(required)) x <- need_columns(x, required[[nm]], nm)
    for (id in intersect(c("patient_id", "hospitalization_id"), names(x))) x[[id]] <- as.character(x[[id]])
    for (ts in grep("_dttm$", names(x), value = TRUE)) x[[ts]] <- utc(x[[ts]])
    if (!is.null(cohort_ids) && "hospitalization_id" %in% names(x)) x <- x[x$hospitalization_id %in% cohort_ids, ]
    else if (!is.null(patient_ids) && "patient_id" %in% names(x)) x <- x[x$patient_id %in% patient_ids, ]
    message("Loaded ", nm, ": ", nrow(x), " rows")
    out[[nm]] <- tibble::as_tibble(x)
  }
  out
}

init_run <- function(cfg) {
  root <- file.path(cfg$repo, "output", paste0("run_", cfg$site_name, "_", cfg$analysis$run_id))
  dir.create(file.path(root, "local"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "share"), recursive = TRUE, showWarnings = FALSE)
  root
}

# Exact aggregate release is the project default. Legacy mode exists solely to
# reproduce historical disclosure-filtered exports; no patient rows are allowed.
aggregate_exact <- function() identical(getOption('clif.aggregate_policy','exact'),'exact')
aggregate_min_cell <- function(requested=10L) if(aggregate_exact())0L else requested
aggregate_count_band <- function(n,width=10L) {
 if(aggregate_exact())return(c(lower=n,upper=n))
 c(lower=width*floor(n/width),upper=width*floor(n/width)+width-1)
}
promote_exact_counts <- function(x) {
 if(!aggregate_exact())return(x)
 for(base in c('n','events','patients','observed','missing','at_risk','classified')){
  lo<-paste0(base,'_lower');hi<-paste0(base,'_upper')
  if(all(c(lo,hi)%in%names(x)) && !base%in%names(x) && all(is.na(x[[lo]]) | is.na(x[[hi]]) | x[[lo]]==x[[hi]]))x[[base]]<-x[[lo]]
 }
 x
}

# Historical mode: each partition is fully masked if any constituent count is small. This avoids
# recovering a suppressed cell from the other cells and a known denominator.
# Central code must not treat suppressed numeric values as zero.
protect_aggregate <- function(x, count_cols, partition, min_cell, value_cols = NULL) {
  if (!nrow(x)) return(x)
  count_cols <- intersect(count_cols, names(x))
  x$cell_status <- "observed"
  if(aggregate_exact())return(x)
  bad <- rep(FALSE, nrow(x))
  for (nm in count_cols) bad <- bad | (!is.na(x[[nm]]) & x[[nm]] > 0 & x[[nm]] < min_cell)
  if (length(partition)) {
    key <- do.call(paste, c(x[partition], sep = "\r"))
    bad <- key %in% unique(key[bad])
  }
  if (is.null(value_cols)) value_cols <- names(x)[vapply(x, is.numeric, logical(1))]
  for (nm in intersect(value_cols, names(x))) x[[nm]][bad] <- NA_real_
  x$cell_status[bad] <- "suppressed_partition"
  x
}

write_share <- function(x, root, name, cfg) {
  x <- promote_exact_counts(x)
  forbidden <- intersect(names(x), c("patient_id", "hospitalization_id", "t0", "event_time", "admission_dttm", "discharge_dttm"))
  if (length(forbidden)) stop("Patient-level columns in aggregate export: ", paste(forbidden, collapse = ", "))
  if (any(vapply(x, inherits, logical(1), "POSIXt"))) stop("Raw timestamps in aggregate export")
  x$schema_version <- rep(RUN_SCHEMA, nrow(x))
  x$site <- rep(cfg$site_name, nrow(x))
  x$run_id <- rep(cfg$analysis$run_id, nrow(x))
  readr::write_csv(x, file.path(root, "share", paste0(name, ".csv")), na = "")
}

run_manifest <- function(cfg, root, phase) {
  files <- list.files("code", recursive = TRUE, full.names = TRUE, pattern = "\\.R$")
  manifest <- list(schema_version = RUN_SCHEMA, site = cfg$site_name, run_id = cfg$analysis$run_id,
    phase = phase, parameters = cfg$analysis, analysis_plan = analysis_specifications(cfg),
    code_md5 = as.list(tools::md5sum(files)),
    R_version = R.version.string,
    package_versions = lapply(stats::setNames(c("dplyr", "TraMineR", "survival", "metafor", "comorbidity"),
      c("dplyr", "TraMineR", "survival", "metafor", "comorbidity")), function(p) as.character(utils::packageVersion(p))))
  jsonlite::write_json(manifest, file.path(root, "share", paste0("manifest_", phase, ".json")), auto_unbox = TRUE, pretty = TRUE, null = "null")
}

code_fingerprint <- function() {
  files <- sort(list.files("code", recursive = TRUE, full.names = TRUE, pattern = "\\.R$"))
  digest::digest(as.list(tools::md5sum(files)), algo = "sha256")
}
