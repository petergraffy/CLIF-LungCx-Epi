association_terms <- function(d) {
  candidates <- c("age_years", "sex_category", "race_category", "admission_year",
    "sofa_total", "charlson_noncancer", "advanced_cancer", "code_at_admission")
  candidates[vapply(d[candidates], function(x) length(unique(x[!is.na(x)])) > 1L, logical(1))]
}

fit_associations <- function(d, reference = "P1") {
  d$phenotype <- stats::relevel(factor(d$phenotype), ref = reference)
  # Use one common complete-case population for all three endpoints.
  baseline <- c("age_years", "sex_category", "race_category", "admission_year", "sofa_total",
                "charlson_noncancer", "advanced_cancer", "code_at_admission")
  keep <- stats::complete.cases(d[c(baseline, "phenotype", "time_from_landmark", "event_type", "patient_id")]) & d$time_from_landmark > 0
  audit <- tibble::tibble(variable = c(baseline, "phenotype"), missing_n = vapply(d[c(baseline, "phenotype")], function(x) sum(is.na(x)), integer(1)))
  x <- droplevels(d[keep, ])
  metadata <- list(); estimates <- list(); covariance <- list()
  for (endpoint in c("death", "hospice", "death_or_hospice")) {
    event <- if (endpoint == "death_or_hospice") x$event_type %in% c("death", "hospice") else x$event_type == endpoint
    x$event <- as.integer(event)
    formula <- stats::as.formula(paste("survival::Surv(time_from_landmark, event) ~ phenotype +",
       paste(c(association_terms(x), "cluster(patient_id)"), collapse = " + ")))
    warnings <- character()
    fit <- if (nrow(x) && length(unique(x$phenotype)) > 1 && sum(event) > 0) tryCatch(
      withCallingHandlers(survival::coxph(formula, data = x, ties = "efron", x = TRUE, model = TRUE),
        warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }), error = identity) else simpleError("Insufficient patients, phenotypes, or events")
    ok <- inherits(fit, "coxph") && !length(warnings) && all(is.finite(stats::coef(fit)))
    metadata[[endpoint]] <- tibble::tibble(endpoint, candidate_n = nrow(d), fitted_n = if (inherits(fit, "coxph")) fit$n else nrow(x),
      patient_n = dplyr::n_distinct(x$patient_id), events = sum(event), non_events = sum(!event), excluded_n = sum(!keep),
      status = if (ok) "ok" else if (inherits(fit, "error")) conditionMessage(fit) else paste(c("nonestimable_or_warning", warnings), collapse = "; "),
      formula = paste(deparse(formula), collapse = " "), estimand = "cause_specific_hazard")
    if (!ok) next
    cf <- stats::coef(fit); vc <- stats::vcov(fit)
    ix <- grep("^phenotype", names(cf))
    zph <- tryCatch(survival::cox.zph(fit)$table, error = function(e) NULL)
    estimates[[endpoint]] <- tibble::tibble(endpoint, term = names(cf)[ix],
      phenotype = sub("^phenotype", "", names(cf)[ix]), reference,
      estimate = unname(cf[ix]), std_error = sqrt(diag(vc)[ix]),
      fitted_n = fit$n, events = fit$nevent,
      ph_global_p = if (!is.null(zph)) zph["GLOBAL", "p"] else NA_real_,
      effect_type = "cause_specific_HR")
    vv <- as.data.frame(as.table(vc[ix, ix, drop = FALSE]))
    names(vv) <- c("term1", "term2", "covariance")
    vv$endpoint <- endpoint; covariance[[endpoint]] <- vv
  }
  list(estimates = dplyr::bind_rows(estimates), metadata = dplyr::bind_rows(metadata),
       covariance = dplyr::bind_rows(covariance), missingness = audit)
}

aj_grid <- function(d, days = c(0, 1, 3, 7, 14, 28)) {
  out <- list()
  for (ph in sort(unique(d$phenotype[!is.na(d$phenotype)]))) {
    z <- d[d$phenotype == ph & !is.na(d$phenotype) & d$time_from_landmark > 0, ]
    z$event <- factor(z$event_type, levels = c("censored", "death", "hospice", "live_discharge"))
    # Each hospitalization restarts its own landmark clock. A patient may have
    # multiple episodes; cluster their uncertainty without merging episodes.
    fit <- survival::survfit(survival::Surv(time_from_landmark, event) ~ 1,
      data = z, id = hospitalization_id, cluster = patient_id, robust = TRUE)
    # Do not extend the curve beyond the last observed follow-up.
    times <- days[days <= max(z$time_from_landmark)]
    sm <- summary(fit, times = times)
    out[[ph]] <- dplyr::bind_rows(lapply(seq_along(fit$states), function(j) {
      tibble::tibble(phenotype = ph, day = sm$time, event = fit$states[j],
        cumulative_incidence = sm$pstate[, j], std_error = sm$std.err[, j],
        n_at_risk = if (is.matrix(sm$n.risk)) sm$n.risk[, 1] else sm$n.risk, denominator = nrow(z))
    }))
  }
  dplyr::bind_rows(out)
}

summarize_clusters <- function(d, cluster = "local_cluster") {
  vars <- intersect(c("age_years", "sofa_total", "sofa_observed_domains", "charlson_noncancer",
    "hours_IMV", "hours_HFNC", "hours_CPAP", "hours_NIPPV", "vasopressor_documented_hours",
    "support_coverage", "physiology_unknown_hours"), names(d))
  continuous <- dplyr::bind_rows(lapply(vars, function(v) {
    d |> dplyr::group_by(.data[[cluster]]) |> dplyr::summarise(variable = v,
      n = sum(!is.na(.data[[v]])), missing_n = sum(is.na(.data[[v]])),
      sum = sum(.data[[v]], na.rm = TRUE), sum_squares = sum(.data[[v]]^2, na.rm = TRUE),
      mean = if (all(is.na(.data[[v]]))) NA_real_ else mean(.data[[v]], na.rm = TRUE),
      .groups = "drop")
  }))
  categories <- intersect(c("death_in_hosp", "hospice_discharge", "death_or_hospice", "any_imv", "any_vasopressor",
    "advanced_cancer", "metastatic_cancer", "malignant_pleural_effusion", "cachexia", "code_at_admission",
    "code_at_landmark", "new_dni", "new_comfort", "sex_category", "race_category", "first_hospitalization",
    "landmark_reason", "histology", "stage", "ecog", "therapy_line", "therapy_intent"), names(d))
  categorical <- dplyr::bind_rows(lapply(categories, function(v) {
    z <- d; z$value <- ifelse(is.na(z[[v]]), "MISSING", as.character(z[[v]]))
    z |> dplyr::group_by(.data[[cluster]]) |> dplyr::mutate(denominator = dplyr::n()) |> dplyr::ungroup() |>
      dplyr::count(.data[[cluster]], denominator, value, name = "n") |>
      dplyr::mutate(variable = v, percent = 100 * n / denominator)
  }))
  list(continuous = continuous, categorical = categorical)
}
