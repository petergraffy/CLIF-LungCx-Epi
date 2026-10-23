prepare_revision <- function(cfg) {
  root <- init_run(cfg)
  cache <- file.path(root, "local", "prepared.rds")
  if (file.exists(cache)) {
    old <- readRDS(cache)
    if (!identical(old$parameters, cfg$analysis) || !identical(old$source_id, code_fingerprint()))
      stop("Cached preparation parameters/source differ; choose a new run_id")
    return(old)
  }
  tables <- load_clif(cfg)
  built <- build_cohort(tables, cfg); d <- built$cohort
  for (nm in names(tables)) {
    if ("hospitalization_id" %in% names(tables[[nm]])) tables[[nm]] <- tables[[nm]] |>
      dplyr::filter(hospitalization_id %in% d$hospitalization_id)
    else if ("patient_id" %in% names(tables[[nm]])) tables[[nm]] <- tables[[nm]] |>
      dplyr::filter(patient_id %in% d$patient_id)
  }
  d <- cancer_features(d, tables$hospital_diagnosis) |>
    dplyr::left_join(sofa_features(d, tables, cfg), by = "hospitalization_id")
  # Optional site extension: standardized variables, kept local and descriptive.
  if (!is.null(tables$oncology)) {
    need_columns(tables$oncology, "hospitalization_id", "oncology")
    if (anyDuplicated(tables$oncology$hospitalization_id)) stop("oncology requires one row per hospitalization")
    allowed <- intersect(c("hospitalization_id", "histology", "stage", "ecog", "therapy_line", "therapy_intent"), names(tables$oncology))
    d <- dplyr::left_join(d, tables$oncology[allowed], by = "hospitalization_id")
  }
  write_share(protect_aggregate(built$flow, "n", character(), cfg$analysis$min_cell), root, "cohort_exclusions", cfg)
  disposition <- d |> dplyr::count(discharge_category, event_type, valid_outcome, event_time_source, name = "n")
  write_share(protect_aggregate(disposition, "n", character(), cfg$analysis$min_cell), root, "outcome_timestamp_audit", cfg)
  available <- dplyr::bind_rows(lapply(c("code_status", "oncology"), function(nm) {
    tibble::tibble(source = nm, available = !is.null(tables[[nm]]), rows = if (is.null(tables[[nm]])) 0L else nrow(tables[[nm]]))
  }))
  write_share(protect_aggregate(available, "rows", "source", cfg$analysis$min_cell), root, "source_availability", cfg)
  annual <- d |> dplyr::group_by(admission_year) |>
    dplyr::summarise(n = dplyr::n(), outcome_known_n = sum(valid_outcome),
      death_n = sum(death_in_hosp, na.rm = TRUE), hospice_n = sum(hospice_discharge, na.rm = TRUE),
      live_discharge_n = sum(event_type == "live_discharge" & valid_outcome), unknown_n = sum(!valid_outcome),
      .groups = "drop")
  write_share(protect_aggregate(annual, c("n", "death_n", "hospice_n", "live_discharge_n", "unknown_n"),
    "admission_year", cfg$analysis$min_cell, setdiff(names(annual), "admission_year")), root, "cohort_by_year", cfg)
  patients <- tibble::tibble(hospitalizations = nrow(d), patients = dplyr::n_distinct(d$patient_id),
    repeat_hospitalizations = sum(!d$first_hospitalization))
  write_share(protect_aggregate(patients, names(patients), character(), cfg$analysis$min_cell), root, "cohort_patient_counts", cfg)
  prepared <- list(cohort = d, tables = tables, parameters = cfg$analysis, root = root, source_id = code_fingerprint())
  saveRDS(prepared, cache)
  run_manifest(cfg, root, "prepare")
  prepared
}

derive_revision <- function(prepared, cfg) {
  root <- prepared$root
  cache <- file.path(root, "local", "derivation.rds")
  if (file.exists(cache)) stop("Derivation already exists; select a new run_id for a fresh derivation")
  all <- prepared$cohort; tables <- prepared$tables
  specifications <- analysis_specifications(cfg)
  write_share(specifications, root, "analysis_plan", cfg)
  onset_detection <- first_detected_arf(all, tables, cfg)
  saveRDS(onset_detection, file.path(root, "local", "arf_detection.rds"))
  detection_counts <- onset_detection |> dplyr::count(arf_detection_status, arf_initial_evidence, name = "n")
  write_share(protect_aggregate(detection_counts, "n", character(), cfg$analysis$min_cell),
    root, "arf_detection_audit", cfg)
  analyses <- list(); issues <- list()
  for (s in seq_len(nrow(specifications))) {
    window <- specifications$window[s]; variant <- specifications$variant[s]
    analysis_id <- specifications$analysis_id[s]
    message("Deriving ", analysis_id)
    d <- if (variant == "arf_onset") arf_onset_cohort(all, onset_detection) else landmark_cohort(all, window)
    flow <- d |> dplyr::count(landmark_reason, name = "n") |> dplyr::mutate(analysis_id = analysis_id)
    write_share(protect_aggregate(flow, "n", "analysis_id", cfg$analysis$min_cell), root, paste0(analysis_id, "_landmark_flow"), cfg)
    if (variant == "absorbing") d <- d[d$valid_outcome, ] else d <- d[d$landmark_reason == "eligible", ]
    if (nrow(d) < 4) {
      issues[[analysis_id]] <- tibble::tibble(analysis_id, status = "insufficient_landmark_eligible_hospitalizations")
      next
    }
    panel <- build_panel(d, tables, window, cfg, absorbing = variant == "absorbing",
      estimate_low_flow = variant == "low_flow_fio2",
      carry = switch(variant, arf_carry_0 = 0, arf_carry_12 = 12, NULL))
    summaries <- support_summaries(panel)
    d <- d |> dplyr::left_join(summaries, by = "hospitalization_id") |>
      dplyr::left_join(code_status_features(d, tables$code_status, window, cfg), by = "hospitalization_id") |>
      dplyr::left_join(vasopressor_features(d, tables$medication_admin_continuous, window), by = "hospitalization_id")
    # Every eligible hospitalization contributes a complete hourly sequence.
    # Missing support remains UNKNOWN, including entirely unobserved sequences.
    # An event in the first hour can yield an entirely terminal sequence.
    d$sequence_reason <- "included"
    qc <- d |> dplyr::count(sequence_reason, name = "n") |> dplyr::mutate(analysis_id = analysis_id)
    write_share(protect_aggregate(qc, "n", "analysis_id", cfg$analysis$min_cell), root, paste0(analysis_id, "_sequence_flow"), cfg)
    field <- if (variant == "support_only") "support" else "joint"
    fit <- tryCatch(sequence_fit(panel, field, cfg$analysis$k_max), error = identity)
    if (inherits(fit, "error")) {
      issues[[analysis_id]] <- tibble::tibble(analysis_id, status = conditionMessage(fit)); next
    }
    assignments <- tibble::tibble(hospitalization_id = fit$ids, local_cluster = fit$labels)
    stopifnot(nrow(assignments) == nrow(d), !anyDuplicated(assignments$hospitalization_id),
      setequal(assignments$hospitalization_id, d$hospitalization_id), !anyNA(assignments$local_cluster))
    d <- d |> dplyr::inner_join(assignments, by = "hospitalization_id")
    sig <- cluster_signature(panel, assignments, field) |> dplyr::mutate(analysis_id = analysis_id)
    sig <- protect_aggregate(sig, "n", c("analysis_id", "local_cluster", "h"), cfg$analysis$min_cell,
                              c("n", "probability"))
    # Small cluster sizes are themselves suppressed, along with all profile values.
    small <- sig$cluster_n < cfg$analysis$min_cell
    for (v in c("n", "cluster_n", "probability")) sig[[v]][small] <- NA_real_
    sig$cell_status[small] <- "suppressed_cluster"
    write_share(sig, root, paste0(analysis_id, "_signature"), cfg)
    write_share(fit$silhouette |> dplyr::mutate(analysis_id = analysis_id), root, paste0(analysis_id, "_silhouette"), cfg)
    summary <- summarize_clusters(d)
    for (kind in names(summary)) {
      z <- summary[[kind]] |> dplyr::mutate(analysis_id = analysis_id)
      counts <- if (kind == "continuous") c("n", "missing_n") else "n"
      z <- protect_aggregate(z, counts, c("analysis_id", "local_cluster", "variable"), cfg$analysis$min_cell,
        setdiff(names(z)[vapply(z, is.numeric, logical(1))], "local_cluster"))
      write_share(z, root, paste0(analysis_id, "_", kind), cfg)
    }
    by_code <- d |> dplyr::group_by(local_cluster, code_at_admission) |>
      dplyr::summarise(n = dplyr::n(), death_n = sum(death_in_hosp), hospice_n = sum(hospice_discharge),
        live_discharge_n = sum(event_type == "live_discharge"),
        composite_n = sum(death_or_hospice), .groups = "drop")
    by_code <- protect_aggregate(by_code, c("n", "death_n", "hospice_n", "live_discharge_n"),
      "local_cluster", cfg$analysis$min_cell, c("n", "death_n", "hospice_n", "live_discharge_n", "composite_n"))
    write_share(by_code, root, paste0(analysis_id, "_outcomes_by_code_status"), cfg)
    missing <- d |> dplyr::group_by(local_cluster) |>
      dplyr::summarise(n = dplyr::n(), support_coverage = mean(support_coverage),
        carried_hours = sum(carried_hours), physiology_unknown_hours = sum(physiology_unknown_hours),
        sofa_missing_n = sum(is.na(sofa_total)), .groups = "drop")
    write_share(protect_aggregate(missing, c("n", "sofa_missing_n"), "local_cluster", cfg$analysis$min_cell,
      setdiff(names(missing), "local_cluster")), root, paste0(analysis_id, "_missingness"), cfg)
    if (variant == "primary" && cfg$analysis$bootstrap_reps > 0) {
      stability <- bootstrap_stability(fit, d$patient_id[match(fit$ids, d$hospitalization_id)], cfg$analysis$bootstrap_reps, cfg$analysis$seed)
      if (nrow(stability)) {
        stability <- stability |> dplyr::group_by(local_cluster) |>
          dplyr::summarise(completed_replicates = dplyr::n(), mean_jaccard = mean(jaccard),
            lower_jaccard = stats::quantile(jaccard, .025), mean_adjusted_rand = mean(adjusted_rand, na.rm = TRUE), .groups = "drop")
        write_share(stability, root, paste0(analysis_id, "_stability"), cfg)
      }
    }
    # Keep sequence matrices, clinical data, and patient IDs strictly local.
    fit$distance <- NULL
    analyses[[analysis_id]] <- list(data = d, fit = fit, window = window, variant = variant,
      specification = specifications[s, ])
    issues[[analysis_id]] <- tibble::tibble(analysis_id, status = "ok")
  }
  statuses <- dplyr::bind_rows(issues)
  write_share(statuses, root, "derivation_status", cfg)
  saveRDS(list(analyses = analyses, parameters = cfg$analysis, source_id = code_fingerprint()), cache)
  run_manifest(cfg, root, "derive")
  if (!"h72_primary" %in% names(analyses)) stop("Primary derivation failed; inspect derivation_status.csv")
  invisible(analyses)
}

model_revision <- function(cfg, mapping_path) {
  root <- init_run(cfg)
  cache <- readRDS(file.path(root, "local", "derivation.rds"))
  if (!identical(cache$parameters, cfg$analysis) || !identical(cache$source_id, code_fingerprint()))
    stop("Model parameters/source differ from derivation")
  mapping <- readr::read_csv(mapping_path, show_col_types = FALSE)
  need_columns(mapping, c("site", "run_id", "schema_version", "analysis_id", "local_cluster", "phenotype", "reference", "consensus_id", "source_id"), "mapping")
  if (any(mapping$schema_version != RUN_SCHEMA) || anyNA(mapping$schema_version)) stop("Incompatible mapping schema")
  mapping <- mapping[mapping$site == cfg$site_name & mapping$run_id == cfg$analysis$run_id, ]
  if (!nrow(mapping)) stop("Mapping does not contain this site and run_id")
  if (any(mapping$source_id != cache$source_id)) stop("Mapping source version differs from local derivation")
  if (anyDuplicated(mapping[c("analysis_id", "local_cluster")])) stop("Duplicate mapping keys")
  if (length(unique(mapping$consensus_id)) != 1) stop("Ambiguous consensus version")
  consensus_id <- unique(mapping$consensus_id)
  model_dir <- file.path(root, "share", paste0("models_", consensus_id))
  if (dir.exists(model_dir)) stop("Model export already exists for this consensus")
  dir.create(model_dir)
  # write_share writes to <root>/share; use an isolated model package root.
  model_root <- file.path(root, "local", paste0("model_package_", consensus_id))
  dir.create(file.path(model_root, "share"), recursive = TRUE)
  for (id in names(cache$analyses)) {
    analysis <- cache$analyses[[id]]
    if (analysis$specification$outcome_modeling == "descriptive_only") next
    map <- mapping[mapping$analysis_id == id, ]
    if (!nrow(map)) next
    d <- analysis$data |> dplyr::left_join(map |> dplyr::select(local_cluster, phenotype), by = "local_cluster")
    reference <- unique(map$reference)
    if (length(reference) != 1) stop("Ambiguous reference for ", id)
    alignment <- d |> dplyr::count(local_cluster, phenotype, name = "n")
    write_share(protect_aggregate(alignment, "n", character(), cfg$analysis$min_cell), model_root, paste0(id, "_mapping_coverage"), cfg)
    subsets <- list(all = rep(TRUE, nrow(d)), first_hospitalization = d$first_hospitalization,
      full_code_at_admission = d$code_at_admission == "full", no_covid = !d$covid,
      arf_or_high_support = d$high_support | d$measured_arf,
      sofa_missing_zero_sensitivity = rep(TRUE, nrow(d)))
    for (subset in names(subsets)) {
      x <- d[subsets[[subset]], ]
      if (subset == "sofa_missing_zero_sensitivity") x$sofa_total <- x$sofa_missing_zero_sensitivity
      if (!reference %in% x$phenotype) next
      assoc <- fit_associations(x, reference)
      for (kind in names(assoc)) {
        z <- assoc[[kind]]
        if (!nrow(z)) next
        z$analysis_id <- id; z$subset <- subset; z$consensus_id <- consensus_id
        z$analysis_role <- analysis$specification$analysis_role
        z$trajectory_origin <- analysis$specification$trajectory_origin
        z$landmark_hours <- analysis$window
        if (kind == "metadata") z <- protect_aggregate(z, c("fitted_n", "events", "non_events"), c("analysis_id", "subset"), cfg$analysis$min_cell,
          intersect(c("candidate_n", "fitted_n", "patient_n", "events", "non_events", "excluded_n"), names(z)))
        if (kind %in% c("estimates", "covariance")) {
          unsafe <- assoc$metadata$endpoint[assoc$metadata$events < cfg$analysis$min_cell | assoc$metadata$non_events < cfg$analysis$min_cell]
          z <- z[!z$endpoint %in% unsafe, ]
        }
        if (kind == "missingness") z <- protect_aggregate(z, "missing_n", c("analysis_id", "subset"), cfg$analysis$min_cell, "missing_n")
        write_share(z, model_root, paste0(id, "_", subset, "_", kind), cfg)
      }
    }
    cif <- aj_grid(d[!is.na(d$phenotype), ])
    # CIFs and at-risk counts can disclose small event cells; do not export
    # curves for a phenotype with a small nonzero terminal-event category.
    risks <- d |> dplyr::filter(!is.na(phenotype)) |> dplyr::count(phenotype, event_type)
    unsafe <- unique(risks$phenotype[risks$n > 0 & risks$n < cfg$analysis$min_cell])
    cif <- cif[!cif$phenotype %in% unsafe, ]
    cif <- protect_aggregate(cif, "n_at_risk", "phenotype", cfg$analysis$min_cell,
      c("cumulative_incidence", "std_error", "n_at_risk", "denominator"))
    write_share(cif, model_root, paste0(id, "_cumulative_incidence"), cfg)
    if (analysis$variant == "primary") {
      seq <- analysis$fit$sequence[match(d$hospitalization_id, analysis$fit$ids), , drop = FALSE]
      for (severity in c("measured_components", "missing_zero_sensitivity")) {
      validation_data <- d
      if (severity == "missing_zero_sensitivity") validation_data$sofa_total <- validation_data$sofa_missing_zero_sensitivity
      validation <- validate_nested(validation_data, seq, cfg)
      for (kind in names(validation)) {
        z <- validation[[kind]]; if (!nrow(z)) next
        z$analysis_id <- id; z$consensus_id <- consensus_id
        if (kind == "calibration") z <- protect_aggregate(z, c("events", "non_events"), "model", cfg$analysis$min_cell,
          c("n", "events", "non_events", "mean_predicted", "observed_risk"))
        if (kind == "decision") z <- protect_aggregate(z, c("tp", "fp"), "model", cfg$analysis$min_cell,
          c("n", "tp", "fp", "net_benefit", "treat_all"))
        if (kind == "audit") z <- protect_aggregate(z, c("events", "non_events"), "analysis_id", cfg$analysis$min_cell,
          intersect(c("candidate_n", "complete_case_n", "excluded_n", "events", "non_events"), names(z)))
        if (kind %in% c("metrics", "differences") && !is.null(validation$audit$events) &&
            min(validation$audit$events, validation$audit$non_events) < cfg$analysis$min_cell) next
        write_share(z, model_root, paste0(id, "_validation_", severity, "_", kind), cfg)
      }
      }
    }
  }
  files <- list.files(file.path(model_root, "share"), full.names = TRUE)
  if (!all(file.copy(files, model_dir))) stop("Could not finalize model package")
  file.copy(mapping_path, file.path(model_dir, "consensus_mapping.csv"))
  run_manifest(cfg, root, "model")
  file.copy(file.path(root, "share", "manifest_model.json"), file.path(model_dir, "manifest_model.json"))
  invisible(model_dir)
}
