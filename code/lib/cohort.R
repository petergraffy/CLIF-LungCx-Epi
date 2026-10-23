build_cohort <- function(tables, cfg) {
  p <- with_columns(tables$patient, c("death_dttm", "ethnicity_category"))
  h <- tables$hospitalization
  stopifnot(!anyDuplicated(p$patient_id), !anyDuplicated(h$hospitalization_id),
            !anyNA(p$patient_id), !anyNA(h$patient_id), !anyNA(h$hospitalization_id))
  p$death_dttm <- utc(p$death_dttm)
  dx <- tables$hospital_diagnosis
  code <- toupper(gsub("[^A-Za-z0-9]", "", dx$diagnosis_code))
  fmt <- toupper(dx$diagnosis_code_format)
  lung <- dx$hospitalization_id[which(num(dx$poa_present) == 1 &
    ((grepl("10", fmt) & grepl("^C34", code)) | (grepl("9", fmt) & grepl("^162", code))))]
  icu <- tables$adt |> dplyr::filter(tolower(location_category) == "icu", !is.na(in_dttm))
  bounds <- icu |> dplyr::group_by(hospitalization_id) |>
    dplyr::summarise(t0 = min(in_dttm), .groups = "drop")
  d <- h |> dplyr::left_join(p |> dplyr::select(patient_id, sex_category, race_category, ethnicity_category, death_dttm), by = "patient_id") |>
    dplyr::left_join(bounds, by = "hospitalization_id") |>
    dplyr::mutate(age_years = num(age_at_admission),
      admission_year = lubridate::year(admission_dttm),
      exclusion = dplyr::case_when(
        is.na(admission_dttm) ~ "missing_admission_time",
        admission_dttm < utc(cfg$analysis$start_date) | admission_dttm >= utc(cfg$analysis$end_date) ~ "outside_study_window",
        is.na(age_years) | age_years < 18 ~ "not_known_adult",
        !hospitalization_id %in% lung ~ "no_primary_lung_cancer_poa",
        is.na(t0) ~ "no_icu_admission",
        t0 < admission_dttm | (!is.na(discharge_dttm) & t0 >= discharge_dttm) ~ "invalid_icu_time",
        TRUE ~ "included"))
  flow <- d |> dplyr::count(exclusion, name = "n")
  d <- d |> dplyr::filter(exclusion == "included")
  if (!nrow(d)) stop("No eligible lung cancer ICU hospitalizations")
  disposition <- tolower(trimws(gsub("_", " ", d$discharge_category)))
  d$event_type <- dplyr::case_when(
    grepl("expired|death|died", disposition) ~ "death",
    grepl("hospice", disposition) ~ "hospice",
    disposition %in% c("home", "home health", "skilled nursing facility", "skilled nursing facility (snf)", "snf",
      "rehabilitation", "rehab", "acute inpatient rehab facility", "long term care", "long term care hospital (ltach)", "ltach",
      "acute care hospital", "other hospital", "alive", "discharged alive", "against medical advice", "against medical advice (ama)", "ama",
      "group home", "chemical dependency", "assisted living", "psychiatric hospital", "shelter", "jail") ~ "live_discharge",
    TRUE ~ "unknown")
  # Do not borrow a later outpatient death or a last-vital timestamp as an event.
  valid_death <- d$event_type == "death" & !is.na(d$death_dttm) &
    d$death_dttm >= d$t0 & !is.na(d$discharge_dttm) & d$death_dttm <= d$discharge_dttm
  d$event_time <- d$discharge_dttm
  d$event_time[valid_death] <- d$death_dttm[valid_death]
  d$event_time_source <- ifelse(valid_death, "validated_patient_death", "hospital_discharge")
  d$event_time_source[d$event_type == "unknown"] <- "unknown_disposition"
  d$valid_outcome <- d$event_type != "unknown" & !is.na(d$event_time) & d$event_time > d$t0
  d$death_in_hosp <- ifelse(d$valid_outcome, as.integer(d$event_type == "death"), NA_integer_)
  d$hospice_discharge <- ifelse(d$valid_outcome, as.integer(d$event_type == "hospice"), NA_integer_)
  d$death_or_hospice <- ifelse(d$valid_outcome, as.integer(d$event_type %in% c("death", "hospice")), NA_integer_)
  d$icu_t0 <- d$t0
  d <- d |> dplyr::group_by(patient_id) |> dplyr::arrange(t0, hospitalization_id, .by_group = TRUE) |>
    dplyr::mutate(first_hospitalization = dplyr::row_number() == 1L) |> dplyr::ungroup()
  list(cohort = d, flow = flow, icu = icu)
}

landmark_cohort <- function(d, window) {
  d$landmark <- d$t0 + window * 3600
  d$landmark_reason <- dplyr::case_when(
    !d$valid_outcome ~ "unknown_or_invalid_outcome",
    d$event_time <= d$landmark & d$event_type == "death" ~ "death_in_window",
    d$event_time <= d$landmark & d$event_type == "hospice" ~ "hospice_in_window",
    d$event_time <= d$landmark ~ "live_discharge_in_window",
    TRUE ~ "eligible")
  d$time_from_landmark <- hours_between(d$event_time, d$landmark) / 24
  d
}

code_status_features <- function(d, codes, window, cfg) {
  result <- d |> dplyr::select(hospitalization_id) |>
    dplyr::mutate(code_at_admission = "unknown", code_at_landmark = "unknown",
                  code_changed = NA, new_dni = NA, new_comfort = NA)
  if (is.null(codes) || !nrow(codes)) return(result)
  need_columns(codes, c("patient_id", "start_dttm", "code_status_category"), "code_status")
  # These are distinct categories: DNR alone does not imply DNI.
  raw <- toupper(gsub("[^A-Za-z]", "", codes$code_status_category))
  codes$status <- dplyr::case_when(
    raw %in% c("FULL", "FULLCODE", "PRESUMEFULL", "PRESUMEDFULL") ~ "full",
    raw %in% c("DNRDNI", "DNARDNI", "DNI", "DONOTINTUBATE") ~ "dni",
    raw %in% c("DNR", "DNAR", "UDNR") ~ "dnr_without_documented_dni",
    raw %in% c("AND", "ALLOWNATURALDEATH", "COMFORTMEASURESONLY", "CMO") ~ "comfort",
    TRUE ~ "unknown")
  for (i in seq_len(nrow(d))) {
    baseline_time <- if ("icu_t0" %in% names(d)) d$icu_t0[i] else d$t0[i]
    z <- codes[which(codes$patient_id == d$patient_id[i] & !is.na(codes$start_dttm)), ]
    z <- z[z$start_dttm >= d$admission_dttm[i] - cfg$analysis$code_status_lookback_hours * 3600 &
             z$start_dttm < d$t0[i] + window * 3600, ]
    z <- z[order(z$start_dttm), ]
    admission <- z$status[z$start_dttm <= baseline_time]
    initial <- if (length(admission)) tail(admission, 1) else "unknown"
    after <- z$status[z$start_dttm > baseline_time]
    result$code_at_admission[i] <- initial
    result$code_at_landmark[i] <- if (nrow(z)) tail(z$status, 1) else "unknown"
    result$code_changed[i] <- if (!nrow(z)) NA else any(c(initial, after)[-1] != head(c(initial, after), -1))
    result$new_dni[i] <- if (!nrow(z)) NA else any(after == "dni") & initial != "dni"
    result$new_comfort[i] <- if (!nrow(z)) NA else any(after == "comfort") & initial != "comfort"
  }
  result
}
