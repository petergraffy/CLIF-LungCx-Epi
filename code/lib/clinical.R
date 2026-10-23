cancer_features <- function(d, dx) {
  dx <- dx |> dplyr::filter(num(poa_present) == 1) |>
    dplyr::mutate(code = toupper(gsub("[^A-Za-z0-9]", "", diagnosis_code)),
      version = dplyr::case_when(grepl("10", diagnosis_code_format) ~ "10", grepl("9", diagnosis_code_format) ~ "9", TRUE ~ NA_character_))
  flags <- dx |> dplyr::group_by(hospitalization_id) |>
    dplyr::summarise(
      metastatic_cancer = any((version == "10" & grepl("^C7[789]", code)) | (version == "9" & grepl("^(196|197|198)", code)), na.rm = TRUE),
      malignant_pleural_effusion = any((version == "10" & grepl("^(J910|C785)", code)) | (version == "9" & code == "51181"), na.rm = TRUE),
      cachexia = any((version == "10" & code == "R64") | (version == "9" & code == "7994"), na.rm = TRUE),
      covid = any(version == "10" & code == "U071", na.rm = TRUE), .groups = "drop")
  domains <- list()
  for (v in c("10", "9")) {
    z <- dx |> dplyr::filter(version == v, !is.na(code), nzchar(code)) |>
      dplyr::select(id = hospitalization_id, code)
    if (nrow(z)) domains[[v]] <- comorbidity::comorbidity(z, id = "id", code = "code",
      map = paste0("charlson_icd", v, "_quan"), assign0 = TRUE)
  }
  weights <- c(mi = 1, chf = 1, pvd = 1, cevd = 1, dementia = 1, cpd = 1,
    rheumd = 1, pud = 1, mld = 1, diab = 1, diabwc = 2, hp = 2, rend = 2,
    canc = 2, msld = 3, metacanc = 6, aids = 6)
  scores <- tibble::tibble(hospitalization_id = character(), charlson_noncancer = numeric())
  if (length(domains)) {
    dom <- dplyr::bind_rows(domains) |> dplyr::group_by(id) |>
      dplyr::summarise(dplyr::across(dplyr::all_of(names(weights)), max), .groups = "drop")
    dom$mld[dom$msld == 1] <- 0; dom$diab[dom$diabwc == 1] <- 0
    dom$canc <- 0; dom$metacanc <- 0
    scores <- tibble::tibble(hospitalization_id = dom$id,
      charlson_noncancer = as.numeric(as.matrix(dom[names(weights)]) %*% weights))
  }
  d <- d |> dplyr::left_join(flags, by = "hospitalization_id") |>
    dplyr::left_join(scores, by = "hospitalization_id")
  # Every included encounter has a qualifying POA diagnosis; absent additional
  # coded conditions mean no recorded condition, not confirmed clinical absence.
  for (nm in c("metastatic_cancer", "malignant_pleural_effusion", "cachexia", "covid")) d[[nm]][is.na(d[[nm]])] <- FALSE
  d$advanced_cancer <- as.integer(d$metastatic_cancer | d$malignant_pleural_effusion | d$cachexia)
  d
}

dose_per_kg_min <- function(dose, unit, weight) {
  unit <- tolower(gsub(" ", "", unit)); dose <- num(dose)
  out <- rep(NA_real_, length(dose))
  factor <- c("mcg/kg/min" = 1, "mcg/kg/hr" = 1/60, "mcg/kg/hour" = 1/60,
              "mg/kg/hr" = 1000/60, "mg/kg/hour" = 1000/60)
  for (u in names(factor)) out[which(unit == u)] <- dose[which(unit == u)] * factor[[u]]
  for (u in c("mcg/min", "mg/hr", "mg/hour")) {
    ix <- which(unit == u & !is.na(weight) & weight > 10 & weight < 500)
    out[ix] <- dose[ix] / weight[ix] * ifelse(u == "mcg/min", 1, 1000/60)
  }
  out
}

# Measured-component SOFA proxy: renal criterion uses creatinine, not urine
# output. Missing domains remain NA; no saturation-derived PaO2 or default weight.
sofa_features <- function(d, tables, cfg) {
  rs_by <- split(tables$respiratory_support, tables$respiratory_support$hospitalization_id)
  vit_by <- split(tables$vitals, tables$vitals$hospitalization_id)
  lab_by <- split(tables$labs, tables$labs$hospitalization_id)
  med_by <- split(tables$medication_admin_continuous, tables$medication_admin_continuous$hospitalization_id)
  gcs_by <- split(tables$patient_assessments, tables$patient_assessments$hospitalization_id)
  out <- vector("list", nrow(d))
  max_known <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  for (i in seq_len(nrow(d))) {
    id <- d$hospitalization_id[i]; start <- d$t0[i]; end <- start + 24 * 3600
    if (d$valid_outcome[i]) end <- min(end, d$event_time[i])
    win <- function(x, ts, empty) {
      if (is.null(x)) x <- empty[0, ]
      x[which(!is.na(x[[ts]]) & x[[ts]] >= start & x[[ts]] < end), ]
    }
    v <- win(vit_by[[id]], "recorded_dttm", tables$vitals)
    l <- win(lab_by[[id]], "lab_result_dttm", tables$labs)
    r <- win(rs_by[[id]], "recorded_dttm", tables$respiratory_support)
    m <- win(med_by[[id]], "admin_dttm", tables$medication_admin_continuous)
    g <- win(gcs_by[[id]], "recorded_dttm", tables$patient_assessments)
    lv <- function(cat, lo, hi) {
      x <- num(l$lab_value_numeric[tolower(l$lab_category) == cat]); x[x < lo | x > hi] <- NA_real_; x
    }
    plt <- lv("platelet_count", 0, 3000); bil <- lv("bilirubin_total", 0, 100); cr <- lv("creatinine", 0, 25)
    gc <- num(g$numerical_value[tolower(g$assessment_category) == "gcs_total"]); gc[gc < 3 | gc > 15] <- NA_real_
    mp <- num(v$vital_value[tolower(v$vital_category) == "map"]); mp[mp < 20 | mp > 250] <- NA_real_
    wt <- num(v$vital_value[tolower(v$vital_category) == "weight_kg"]); wt <- wt[!is.na(wt) & wt > 10 & wt < 500]
    weight <- if (length(wt)) stats::median(wt) else NA_real_
    dose <- dose_per_kg_min(m$med_dose, m$med_dose_unit, rep(weight, nrow(m)))
    med <- tolower(m$med_category)
    dop <- max_known(dose[med == "dopamine"]); epi <- max_known(dose[med == "epinephrine"]); norepi <- max_known(dose[med == "norepinephrine"])
    cv <- c(ifelse(mp < 70, 1, 0),
      ifelse(dop > 15, 4, ifelse(dop > 5, 3, ifelse(dop > 0, 2, 0))),
      ifelse(epi > .1, 4, ifelse(epi > 0, 3, 0)), ifelse(norepi > .1, 4, ifelse(norepi > 0, 3, 0)),
      ifelse(dose[med == "dobutamine"] > 0, 2, 0))
    unconverted <- any(num(m$med_dose) > 0 & med %in% c("dopamine", "epinephrine", "norepinephrine", "dobutamine") & is.na(dose), na.rm = TRUE)
    # An unconvertible vasoactive dose prevents a reliable CV score.
    cv_score <- if (unconverted) NA_real_ else max_known(cv)
    po <- physiology_observations(r, v, l, cfg)
    pf <- po$pf; vent <- po$ventilated
    respiratory <- ifelse(pf < 100 & vent, 4, ifelse(pf < 200 & vent, 3,
      ifelse(pf < 300, 2, ifelse(pf < 400, 1, 0))))
    z <- tibble::tibble(hospitalization_id = id,
      sofa_resp = max_known(respiratory), sofa_cv = cv_score,
      sofa_coag = max_known(ifelse(plt < 20, 4, ifelse(plt < 50, 3, ifelse(plt < 100, 2, ifelse(plt < 150, 1, 0))))),
      sofa_liver = max_known(ifelse(bil >= 12, 4, ifelse(bil >= 6, 3, ifelse(bil >= 2, 2, ifelse(bil >= 1.2, 1, 0))))),
      sofa_renal = max_known(ifelse(cr >= 5, 4, ifelse(cr >= 3.5, 3, ifelse(cr >= 2, 2, ifelse(cr >= 1.2, 1, 0))))),
      sofa_cns = max_known(ifelse(gc < 6, 4, ifelse(gc < 10, 3, ifelse(gc < 13, 2, ifelse(gc < 15, 1, 0))))),
      unconverted_cv_dose = unconverted)
    domain <- as.numeric(z[1, grep("^sofa_", names(z))])
    z$sofa_observed_domains <- sum(!is.na(domain))
    z$sofa_total <- if (all(!is.na(domain))) sum(domain) else NA_real_
    z$sofa_missing_zero_sensitivity <- if (all(is.na(domain))) NA_real_ else sum(domain, na.rm = TRUE)
    out[[i]] <- z
  }
  dplyr::bind_rows(out)
}

vasopressor_features <- function(d, meds, window) {
  meds <- meds |> dplyr::filter(tolower(med_category) %in%
    c("norepinephrine", "epinephrine", "dopamine", "phenylephrine", "vasopressin", "angiotensin", "angiotensin_ii"), num(med_dose) > 0) |>
    dplyr::inner_join(d |> dplyr::select(hospitalization_id, t0, event_time), by = "hospitalization_id") |>
    dplyr::filter(admin_dttm >= t0, admin_dttm < t0 + window * 3600,
                  is.na(event_time) | admin_dttm < event_time)
  z <- meds |> dplyr::mutate(h = floor(hours_between(admin_dttm, t0))) |>
    dplyr::group_by(hospitalization_id) |>
    dplyr::summarise(vasopressor_documented_hours = dplyr::n_distinct(h), .groups = "drop")
  d |> dplyr::select(hospitalization_id) |> dplyr::left_join(z, by = "hospitalization_id") |>
    dplyr::mutate(vasopressor_documented_hours = tidyr::replace_na(vasopressor_documented_hours, 0L),
      any_vasopressor = as.integer(vasopressor_documented_hours > 0))
}
