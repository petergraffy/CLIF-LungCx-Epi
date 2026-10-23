# First detected qualifying evidence, not inferred biological onset.
first_detected_arf <- function(d, tables, cfg) {
  rs <- split(tables$respiratory_support, tables$respiratory_support$hospitalization_id)
  vit <- split(tables$vitals, tables$vitals$hospitalization_id)
  lab <- split(tables$labs, tables$labs$hospitalization_id)
  out <- vector("list", nrow(d))
  for (i in seq_len(nrow(d))) {
    id <- d$hospitalization_id[i]; icu <- d$t0[i]
    start <- max(d$admission_dttm[i], icu - cfg$analysis$physiology_carry_hours * 3600)
    end <- icu + 72 * 3600
    if (d$valid_outcome[i]) end <- min(end, d$event_time[i])
    win <- function(x, ts, empty) {
      if (is.null(x)) x <- empty[0, ]
      x[which(!is.na(x[[ts]]) & x[[ts]] >= start & x[[ts]] < end), ]
    }
    obs <- physiology_observations(win(rs[[id]], "recorded_dttm", tables$respiratory_support),
      win(vit[[id]], "recorded_dttm", tables$vitals), win(lab[[id]], "lab_result_dttm", tables$labs), cfg)
    positive <- obs[which(obs$positive & !is.na(obs$available_time) & obs$available_time < end), ]
    first <- if (nrow(positive)) min(positive$available_time) else utc(NA_character_)
    anchor <- if (!is.na(first)) max(icu, first) else utc(NA_character_)
    relation <- if (!nrow(obs)) "no_evaluable_measurements" else if (!nrow(positive)) "no_positive_observed" else
      if (first <= icu) "recent_evidence_before_or_at_icu" else "first_detected_after_icu"
    initial <- if (nrow(positive)) paste(sort(unique(positive$type[positive$available_time == first])), collapse = "_and_") else "unknown"
    out[[i]] <- tibble::tibble(hospitalization_id = id, first_arf_evidence_dttm = first,
      arf_trajectory_origin = anchor, arf_detection_status = relation, arf_initial_evidence = initial,
      arf_detection_delay_hours = hours_between(anchor, icu))
  }
  dplyr::bind_rows(out)
}

arf_onset_cohort <- function(d, detection) {
  d <- dplyr::left_join(d, detection, by = "hospitalization_id")
  d$icu_t0 <- d$t0
  # The sensitivity starts at ICU arrival for recent pre-ICU evidence; it never
  # starts a phenotype window before ICU admission or uses an imputed onset.
  d$t0 <- d$arf_trajectory_origin
  d$landmark <- d$t0 + 72 * 3600
  d$landmark_reason <- dplyr::case_when(
    is.na(d$t0) ~ "no_detected_arf_anchor",
    !d$valid_outcome ~ "unknown_or_invalid_outcome",
    d$event_time <= d$landmark & d$event_type == "death" ~ "death_in_window",
    d$event_time <= d$landmark & d$event_type == "hospice" ~ "hospice_in_window",
    d$event_time <= d$landmark ~ "live_discharge_in_window",
    TRUE ~ "eligible")
  d$time_from_landmark <- hours_between(d$event_time, d$landmark) / 24
  d
}
