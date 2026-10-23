read_return_csv <- function(path) {
  z <- readr::read_csv(path, show_col_types = FALSE)
  need_columns(z, c("schema_version", "site", "run_id"), basename(path))
  if (anyNA(z$schema_version) || any(z$schema_version != RUN_SCHEMA)) stop("Legacy or incompatible return: ", path)
  if (length(unique(z$site)) > 1 || length(unique(z$run_id)) > 1) stop("Mixed sites or runs in ", path)
  z
}

validate_derivation_packages <- function(dirs) {
  manifests <- lapply(dirs, function(dir) {
    path <- file.path(dir, "manifest_derive.json")
    if (!file.exists(path)) stop("Missing derivation manifest: ", dir)
    jsonlite::fromJSON(path, simplifyVector = TRUE)
  })
  if (any(vapply(manifests, function(m) !identical(m$schema_version, RUN_SCHEMA), logical(1)))) stop("Incompatible return schema")
  sites <- vapply(manifests, `[[`, character(1), "site")
  if (anyDuplicated(sites)) stop("More than one return supplied for a site")
  for (m in manifests[-1]) {
    if (!identical(m$parameters, manifests[[1]]$parameters)) stop("Derivation parameters differ between sites")
    if (!identical(m$code_md5, manifests[[1]]$code_md5)) stop("Sites ran different source versions")
  }
  manifests
}

make_consensus_package <- function(dirs, out) {
  if (dir.exists(out)) stop("Consensus destination already exists; use a new directory")
  manifests <- validate_derivation_packages(dirs)
  parameters <- manifests[[1]]$parameters
  plan <- tibble::as_tibble(manifests[[1]]$analysis_plan)
  signatures <- dplyr::bind_rows(lapply(dirs, function(dir) {
    files <- list.files(dir, pattern = "_signature\\.csv$", full.names = TRUE)
    dplyr::bind_rows(lapply(files, read_return_csv))
  }))
  if (!nrow(signatures)) stop("No revision signatures found")
  if (anyDuplicated(signatures[c("site", "analysis_id", "local_cluster", "h", "state")])) stop("Duplicate signature keys")
  fits <- list(); maps <- list(); loso <- list(); issues <- list()
  for (id in sort(unique(signatures$analysis_id))) {
    sig <- signatures[signatures$analysis_id == id, ]
    fit <- tryCatch(fit_consensus(sig, parameters$consensus_k, parameters$consensus_min_overlap), error = identity)
    if (inherits(fit, "error")) {
      issues[[id]] <- tibble::tibble(analysis_id = id, status = conditionMessage(fit)); next
    }
    fits[[id]] <- fit
    maps[[id]] <- fit$mapping |> dplyr::mutate(analysis_id = id, reference = fit$reference,
      matching_inputs = ifelse(grepl("absorbing", id), "joint_respiratory_and_terminal_states_descriptive_only", "hour_specific_respiratory_states_only"))
    loso[[id]] <- consensus_loso(sig, fit) |> dplyr::mutate(analysis_id = id)
    issues[[id]] <- tibble::tibble(analysis_id = id, status = "ok")
  }
  if (!"h72_primary" %in% names(fits)) stop("Cannot derive primary consensus: ", issues$h72_primary$status)
  dir.create(out, recursive = TRUE)
  consensus_id <- substr(digest::digest(list(maps, parameters, manifests[[1]]$code_md5), algo = "sha256"), 1, 16)
  mapping <- dplyr::bind_rows(maps) |>
    dplyr::left_join(plan, by = "analysis_id") |>
    dplyr::left_join(signatures |> dplyr::distinct(site, run_id), by = "site") |>
    dplyr::mutate(schema_version = RUN_SCHEMA, consensus_id = consensus_id,
      source_id = digest::digest(manifests[[1]]$code_md5, algo = "sha256"))
  readr::write_csv(mapping, file.path(out, "consensus_mapping.csv"), na = "")
  readr::write_csv(plan, file.path(out, "analysis_plan.csv"), na = "")
  readr::write_csv(dplyr::bind_rows(loso), file.path(out, "consensus_loso.csv"), na = "")
  readr::write_csv(dplyr::bind_rows(issues), file.path(out, "consensus_status.csv"))
  readr::write_csv(signatures, file.path(out, "source_signatures.csv"), na = "")
  prototypes <- dplyr::bind_rows(lapply(names(fits), function(id) {
    p <- as.data.frame(fits[[id]]$prototypes); p$phenotype <- rownames(p)
    tidyr::pivot_longer(p, -phenotype, names_to = "feature", values_to = "probability") |>
      dplyr::mutate(analysis_id = id)
  }))
  readr::write_csv(prototypes, file.path(out, "consensus_prototypes.csv"), na = "")
  pool_descriptives(dirs, mapping, out)
  saveRDS(fits, file.path(out, "consensus_fits.rds"))
  jsonlite::write_json(list(schema_version = RUN_SCHEMA, consensus_id = consensus_id,
    parameters = parameters, code_md5 = manifests[[1]]$code_md5,
    sites = vapply(manifests, `[[`, character(1), "site")), file.path(out, "consensus_manifest.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  invisible(mapping)
}

pool_descriptives <- function(dirs, mapping, out) {
  strict_sum <- function(x) if (anyNA(x)) NA_real_ else sum(x)
  expected <- mapping |> dplyr::filter(!is.na(phenotype)) |>
    dplyr::group_by(analysis_id, phenotype) |>
    dplyr::summarise(expected_sites = dplyr::n_distinct(site), .groups = "drop")
  keys <- mapping |> dplyr::select(site, run_id, analysis_id, local_cluster, phenotype)
  for (kind in c("continuous", "categorical")) {
    files <- unlist(lapply(dirs, list.files, pattern = paste0("_", kind, "\\.csv$"), full.names = TRUE))
    if (!length(files)) next
    x <- dplyr::bind_rows(lapply(files, read_return_csv)) |>
      dplyr::left_join(keys, by = c("site", "run_id", "analysis_id", "local_cluster")) |>
      dplyr::filter(!is.na(phenotype))
    if (kind == "continuous") {
      result <- x |> dplyr::group_by(analysis_id, phenotype, variable) |>
        dplyr::summarise(contributing_sites = dplyr::n_distinct(site),
          n = strict_sum(n), missing_n = strict_sum(missing_n),
          sum = strict_sum(sum), sum_squares = strict_sum(sum_squares),
          all_observed = all(cell_status == "observed"), .groups = "drop") |>
        dplyr::left_join(expected, by = c("analysis_id", "phenotype")) |>
        dplyr::mutate(complete = all_observed & contributing_sites == expected_sites,
          mean = dplyr::if_else(complete & n > 0, sum / n, NA_real_),
          sd = dplyr::if_else(complete & n > 1, sqrt(pmax(0, (sum_squares - sum^2 / n) / (n - 1))), NA_real_))
      for (nm in c("n", "missing_n", "sum", "sum_squares")) result[[nm]][!result$complete] <- NA_real_
    } else {
      # The denominator is counted once per site/phenotype/variable, rather
      # than once per category. An absent category is zero only when the full
      # site partition is known; suppressed partitions invalidate the total.
      partitions <- x |> dplyr::group_by(site, analysis_id, phenotype, variable) |>
        dplyr::summarise(denominator = if (anyNA(denominator)) NA_real_ else dplyr::first(denominator),
          all_observed = all(cell_status == "observed") && !anyNA(n), .groups = "drop")
      totals <- partitions |> dplyr::group_by(analysis_id, phenotype, variable) |>
        dplyr::summarise(contributing_sites = dplyr::n_distinct(site), denominator = strict_sum(denominator),
          all_observed = all(all_observed), .groups = "drop") |>
        dplyr::left_join(expected, by = c("analysis_id", "phenotype")) |>
        dplyr::mutate(complete = all_observed & contributing_sites == expected_sites)
      result <- x |> dplyr::group_by(analysis_id, phenotype, variable, value) |>
        dplyr::summarise(n = strict_sum(n), .groups = "drop") |>
        dplyr::left_join(totals, by = c("analysis_id", "phenotype", "variable")) |>
        dplyr::mutate(percent = dplyr::if_else(complete & denominator > 0, 100 * n / denominator, NA_real_))
      result$n[!result$complete] <- NA_real_; result$denominator[!result$complete] <- NA_real_
    }
    result$status <- ifelse(result$complete, "complete", "incomplete_suppression_or_missing_site")
    readr::write_csv(result, file.path(out, paste0("pooled_phenotype_", kind, ".csv")), na = "")
  }
  for (stub in c("cohort_by_year", "cohort_patient_counts", "cohort_exclusions", "outcome_timestamp_audit", "source_availability")) {
    files <- file.path(dirs, paste0(stub, ".csv")); files <- files[file.exists(files)]
    if (length(files)) readr::write_csv(dplyr::bind_rows(lapply(files, read_return_csv)),
      file.path(out, paste0("sites_", stub, ".csv")), na = "")
  }
  for (suffix in c("_landmark_flow", "_sequence_flow", "_silhouette", "_stability", "_outcomes_by_code_status")) {
    files <- unlist(lapply(dirs, list.files, pattern = paste0(suffix, "\\.csv$"), full.names = TRUE))
    if (length(files)) readr::write_csv(dplyr::bind_rows(lapply(files, read_return_csv)),
      file.path(out, paste0("sites", suffix, ".csv")), na = "")
  }
  invisible(NULL)
}

pool_effect <- function(d) {
  d <- d[is.finite(d$estimate) & is.finite(d$std_error) & d$std_error > 0, ]
  if (anyDuplicated(d$site)) stop("Duplicate site effects in meta-analysis")
  if (nrow(d) < 2) return(tibble::tibble(k = nrow(d), status = "insufficient_sites"))
  fit <- tryCatch(metafor::rma.uni(yi = estimate, sei = std_error, data = d,
                                  method = "REML", test = "adhoc"), error = identity)
  if (inherits(fit, "error")) return(tibble::tibble(k = nrow(d), status = conditionMessage(fit)))
  pr <- stats::predict(fit)
  tibble::tibble(k = fit$k, status = "ok", estimate = as.numeric(fit$b), std_error = fit$se,
    conf_low = fit$ci.lb, conf_high = fit$ci.ub, tau2 = fit$tau2, i2_percent = fit$I2,
    prediction_low = pr$pi.lb, prediction_high = pr$pi.ub,
    hazard_ratio = exp(as.numeric(fit$b)), hr_low = exp(fit$ci.lb), hr_high = exp(fit$ci.ub),
    hr_prediction_low = exp(pr$pi.lb), hr_prediction_high = exp(pr$pi.ub),
    method = "REML_modified_Hartung_Knapp")
}

pool_model_packages <- function(dirs, out) {
  if (dir.exists(out)) stop("Pooling destination already exists")
  manifests <- lapply(dirs, function(dir) jsonlite::fromJSON(file.path(dir, "manifest_model.json")))
  if (anyDuplicated(vapply(manifests, `[[`, character(1), "site"))) stop("Duplicate site model packages")
  for (m in manifests) if (!identical(m$schema_version, RUN_SCHEMA) ||
    !identical(m$parameters, manifests[[1]]$parameters) || !identical(m$code_md5, manifests[[1]]$code_md5))
      stop("Model packages use incompatible schemas, parameters, or source versions")
  maps <- lapply(dirs, function(dir) {
    x <- readr::read_csv(file.path(dir, "consensus_mapping.csv"), show_col_types = FALSE)
    need_columns(x, c("schema_version", "consensus_id", "source_id"), "mapping")
    if (anyNA(x$schema_version) || any(x$schema_version != RUN_SCHEMA)) stop("Incompatible mapping schema")
    x
  })
  versions <- unique(unlist(lapply(maps, function(z) unique(z$consensus_id))))
  if (length(versions) != 1) stop("Model packages use different consensus versions")
  all <- dplyr::bind_rows(lapply(dirs, function(dir) {
    files <- list.files(dir, pattern = "_estimates\\.csv$", full.names = TRUE)
    dplyr::bind_rows(lapply(files, read_return_csv))
  }))
  if (!nrow(all)) stop("No estimable phenotype effects returned")
  if (any(all$consensus_id != versions)) stop("Incompatible model estimates")
  if (anyDuplicated(all[c("site", "analysis_id", "subset", "endpoint", "phenotype")])) stop("Duplicate site model returns")
  if (any(all$reference != "P1")) stop("Model contrasts have inconsistent reference phenotypes")
  pooled <- all |> dplyr::group_by(analysis_id, subset, endpoint, phenotype, reference, effect_type) |>
    dplyr::group_modify(~ pool_effect(.x)) |> dplyr::ungroup()
  dir.create(out, recursive = TRUE)
  readr::write_csv(all, file.path(out, "site_phenotype_effects.csv"), na = "")
  readr::write_csv(pooled, file.path(out, "pooled_phenotype_effects.csv"), na = "")
  for (pattern in c("_metadata\\.csv$", "_validation_.*\\.csv$", "_cumulative_incidence\\.csv$", "_mapping_coverage\\.csv$")) {
    files <- unlist(lapply(dirs, list.files, pattern = pattern, full.names = TRUE))
    if (!length(files)) next
    # Each export family remains separate: no fabricated denominators, count
    # rescaling, replacement of suppressed cells, or averaging of site medians.
    groups <- split(files, basename(files))
    for (name in names(groups)) {
      x <- dplyr::bind_rows(lapply(groups[[name]], read_return_csv))
      readr::write_csv(x, file.path(out, paste0("sites_", name)), na = "")
    }
  }
  invisible(pooled)
}
