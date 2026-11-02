# Run from repository root with Rscript tests/run_tests.R. Synthetic data only.
source("code/lib/load.R")
source("code/lib/central.R")
check <- function(ok, label) { if (!isTRUE(ok)) stop("FAIL: ", label); cat("PASS:", label, "\n") }
expect_error <- function(expr) inherits(tryCatch({force(expr); NULL}, error = identity), "error")

fixture <- function(n = 160L) {
  set.seed(41)
  t0 <- utc("2020-01-03 00:00:00") + seq_len(n) * 86400
  h <- tibble::tibble(patient_id = paste0("p", seq_len(n)), hospitalization_id = paste0("h", seq_len(n)),
    admission_dttm = t0 - 48 * 3600, discharge_dttm = t0 + sample(96:480, n, TRUE) * 3600,
    age_at_admission = as.character(sample(35:85, n, TRUE)),
    discharge_category = sample(c("Expired", "Hospice", "Home"), n, TRUE, c(.3, .15, .55)))
  p <- tibble::tibble(patient_id = h$patient_id, sex_category = sample(c("Male", "Female"), n, TRUE),
    race_category = sample(c("Black", "White"), n, TRUE), death_dttm = h$discharge_dttm,
    ethnicity_category = "Not Hispanic")
  # Live and hospice cases may have a later outpatient death. It is not a hospital event.
  p$death_dttm[h$discharge_category != "Expired"] <- h$discharge_dttm[h$discharge_category != "Expired"] + 86400 * 10
  adt <- tibble::tibble(hospitalization_id = h$hospitalization_id, in_dttm = t0, out_dttm = t0 + 36 * 3600, location_category = "icu")
  dx <- tibble::tibble(hospitalization_id = h$hospitalization_id, diagnosis_code = "C34.90", diagnosis_code_format = "ICD10CM", poa_present = 1L)
  group <- rep(seq_len(4), length.out = n)
  rs <- dplyr::bind_rows(lapply(seq_len(n), function(i) {
    hh <- 0:71
    state <- switch(as.character(group[i]), "1" = rep("Room Air", 72),
      "2" = ifelse(hh < 24, "IMV", "Nasal Cannula"), "3" = rep("High Flow NC", 72), "4" = rep("IMV", 72))
    tibble::tibble(hospitalization_id = h$hospitalization_id[i], recorded_dttm = t0[i] + hh * 3600 + 60,
      device_category = state, fio2_set = ifelse(state == "Room Air", .21, .4), flow_rate = 4)
  }))
  vit <- dplyr::bind_rows(lapply(seq_len(n), function(i) tibble::tibble(hospitalization_id = h$hospitalization_id[i],
    recorded_dttm = rep(t0[i] + 120, 3), vital_category = c("map", "spo2", "weight_kg"), vital_value = c(80, 98, 80))))
  labs <- dplyr::bind_rows(lapply(seq_len(n), function(i) tibble::tibble(hospitalization_id = h$hospitalization_id[i],
    lab_result_dttm = t0[i] + 120, lab_category = c("po2_arterial", "pco2_arterial", "ph_arterial", "platelet_count", "bilirubin_total", "creatinine"),
    lab_value_numeric = c(90, 40, 7.4, 250, .7, sample(c(.8, 1.5, 2.5), 1)))))
  med <- tibble::tibble(hospitalization_id = character(), admin_dttm = utc(character()), med_category = character(), med_dose = numeric(), med_dose_unit = character())
  gcs <- tibble::tibble(hospitalization_id = h$hospitalization_id, recorded_dttm = t0 + 120, assessment_category = "gcs_total", numerical_value = 15)
  codes <- tibble::tibble(patient_id = h$patient_id, start_dttm = t0 - 3600, code_status_category = "Full")
  list(patient = p, hospitalization = h, adt = adt, hospital_diagnosis = dx, respiratory_support = rs,
       vitals = vit, labs = labs, medication_admin_continuous = med, patient_assessments = gcs, code_status = codes)
}

cfg <- list(site_name = "SyntheticA", repo = getwd(), analysis = run_defaults())
cfg$analysis$bootstrap_reps <- 3L; cfg$analysis$validation_bootstrap_reps <- 10L
cfg$analysis$validation_folds <- 3L; cfg$analysis$k_max <- 4L; cfg$analysis$min_cell <- 0L
tables <- fixture()
built <- build_cohort(tables, cfg); d <- built$cohort
plan <- analysis_specifications(cfg)
check(nrow(plan) == 9L && identical(plan$analysis_id[plan$analysis_role == "primary"], "h72_primary"), "only joint ICU 72-hour analysis is primary")
onset_d <- d[match(paste0("h", 1:4), d$hospitalization_id), ]
det <- first_detected_arf(onset_d, tables, cfg)
check(det$arf_detection_status[1] == "no_positive_observed" && all(det$arf_detection_status[2:4] == "first_detected_after_icu"), "ARF anchor requires observed positive physiology")
anchored <- arf_onset_cohort(onset_d, det)
check(is.na(anchored$t0[1]) && all(anchored$landmark[2:4] == det$arf_trajectory_origin[2:4] + 72 * 3600), "ARF sensitivity has its own origin and 72-hour landmark")
shifted <- onset_d[2, ]; shifted$event_time <- shifted$t0 + 73 * 3600; shifted$event_type <- "death"
shifted_det <- det[2, ]; shifted_det$arf_trajectory_origin <- shifted$t0 + 2 * 3600
check(landmark_cohort(shifted, 72)$landmark_reason == "eligible" &&
  arf_onset_cohort(shifted, shifted_det)$landmark_reason == "death_in_window", "ARF sensitivity excludes events during its shifted trajectory window")
empty <- tables; empty$labs <- empty$labs[0, ]; empty$vitals <- empty$vitals[0, ]
check(first_detected_arf(onset_d[1, ], empty, cfg)$arf_detection_status == "no_evaluable_measurements", "unmeasured ARF is not absence of ARF")
paired <- tables
id2 <- onset_d$hospitalization_id[2]; origin2 <- onset_d$t0[2]
paired$respiratory_support <- paired$respiratory_support |> dplyr::filter(hospitalization_id == id2) |> dplyr::slice(1)
paired$respiratory_support$recorded_dttm <- origin2 + 1800
check(first_detected_arf(onset_d[2, ], paired, cfg)$arf_trajectory_origin == origin2 + 1800, "ARF anchor waits until both paired measurements are available")
pre <- tables
for (nm in c("respiratory_support", "vitals", "labs")) {
  ts <- if (nm == "labs") "lab_result_dttm" else "recorded_dttm"
  use <- pre[[nm]]$hospitalization_id == id2
  pre[[nm]][[ts]][use] <- pre[[nm]][[ts]][use] - 3600
}
check(first_detected_arf(onset_d[2, ], pre, cfg)$arf_trajectory_origin == origin2, "recent pre-ICU ARF evidence anchors at ICU arrival")
later_code <- tables$code_status[2, ]; later_code$start_dttm <- origin2 + 3600; later_code$code_status_category <- "DNR/DNI"
later_d <- onset_d[2, ]; later_d$t0 <- origin2 + 7200
check(code_status_features(later_d, later_code, 72, cfg)$code_at_admission == "unknown", "shifted ARF origin does not redefine ICU admission code status")
check(nrow(d) == nrow(tables$hospitalization), "cohort requires no geography, pollution, or complete demographics")
check(all(d$event_time[d$event_type == "hospice"] == d$discharge_dttm[d$event_type == "hospice"]), "hospice uses discharge time, never later death")
check(all(d$event_time[d$event_type == "live_discharge"] == d$discharge_dttm[d$event_type == "live_discharge"]), "later death cannot become live-discharge event time")
land <- landmark_cohort(d, 72)
check(all(land$landmark == d$t0 + 72 * 3600), "landmark anchored to ICU, after ward transfer")

edge <- d[1:4, ]; edge$event_type <- c("death", "hospice", "live_discharge", "death")
edge$event_time <- edge$t0 + c(6, 24, 48, 72) * 3600
check(all(landmark_cohort(edge, 72)$landmark_reason != "eligible"), "events at or before landmark excluded")
absorb <- build_panel(edge, tables, 72, cfg, absorbing = TRUE)
check(all(absorb$support[absorb$hospitalization_id == edge$hospitalization_id[1] & absorb$h >= 6] == "DEATH"), "death is absorbing; support never extended after death")
check(all(absorb$support[absorb$hospitalization_id == edge$hospitalization_id[2] & absorb$h >= 24] == "HOSPICE"), "hospice absorbing state")
check(all(absorb$support[absorb$hospitalization_id == edge$hospitalization_id[3] & absorb$h >= 48] == "LIVE_DISCHARGE"), "live discharge absorbing state")

panel <- build_panel(d, tables, 72, cfg)
check(nrow(panel) == 72 * nrow(d) && max(panel$h) == 71, "exactly 72 half-open hourly intervals")
check(all(support_state(c("HFNC", "CPAP", "NIPPV")) == c("HFNC", "CPAP", "NIPPV")), "HFNC CPAP and NIPPV remain separate")
missing <- tables; missing$labs <- missing$labs[0, ]; missing$vitals <- missing$vitals[0, ]
mp <- build_panel(d[1, ], missing, 72, cfg)
check(all(mp$physiology == "HU_CU"), "missing physiology never encoded as normal")
late <- tables; late$respiratory_support <- late$respiratory_support |>
  dplyr::filter(hospitalization_id == d$hospitalization_id[1]) |>
  dplyr::slice(11)
lp <- build_panel(d[1, ], late, 72, cfg)
check(all(lp$support[lp$h < 10] == "UNKNOWN") && all(lp$support[lp$h > 17] == "UNKNOWN"), "no backfill and bounded support carry")
none <- tables; none$respiratory_support <- none$respiratory_support[0, ]
np <- build_panel(d[1, ], none, 72, cfg)
check(nrow(np) == 72 && all(np$support == "UNKNOWN"), "patients without support records retained as unknown")
zero_hour <- edge[1, ]; zero_hour$event_time <- zero_hour$t0 + 1800
zero_hour_panel <- build_panel(zero_hour, none, 72, cfg, absorbing = TRUE)
check(nrow(zero_hour_panel) == 72 && all(zero_hour_panel$joint == "DEATH"), "first-hour death retained as a complete absorbing sequence without pre-event bins")
check(nearest_value(utc("2020-01-01 00:00:00"), utc("2020-01-01 02:00:00"), 1, 1) |> is.na(), "FiO2 pairing obeys time tolerance")
check(identical(standard_fio2(c(40, .4, 2, 101)), c(.4, .4, NA_real_, NA_real_)), "FiO2 units and bounds normalized")
check(identical(standard_fio2(c(NA, 40, .4)), c(NA_real_, .4, .4)), "missing FiO2 values are handled without invalid subscripts")

code <- tables$code_status[match(d$patient_id[1:2], tables$code_status$patient_id), ]; code$code_status_category <- c("DNR", "DNR/DNI")
cs <- code_status_features(d[1:2, ], code, 72, cfg)
check(identical(cs$code_at_admission, c("dnr_without_documented_dni", "dni")), "DNR alone never treated as DNI")
stale <- code; stale$start_dttm <- d$admission_dttm[1:2] - 86400
check(all(code_status_features(d[1:2, ], stale, 72, cfg)$code_at_admission == "unknown"), "stale prior-admission code status not propagated")
check(all(is.na(sofa_features(d[1:2, ], missing, cfg)$sofa_total)), "missing SOFA components do not become zero")
cf <- cancer_features(d, tables$hospital_diagnosis)
check(all(cf$charlson_noncancer == 0), "lung cancer does not inflate noncancer Charlson")
sf <- sofa_features(d, tables, cfg)
check(all(sf$sofa_observed_domains == 6) && all(!is.na(sf$sofa_total)), "SOFA proxy uses six measured domains in complete fixture")

fit <- sequence_fit(panel, k_max = 4)
check(length(fit$labels) == nrow(d), "all eligible sequences receive a local cluster")
unknown_panel <- panel
unknown_panel$joint[unknown_panel$hospitalization_id == d$hospitalization_id[1]] <- "UNKNOWN|HU_CU"
unknown_fit <- sequence_fit(unknown_panel, k_max = 4)
check(d$hospitalization_id[1] %in% unknown_fit$ids && length(unknown_fit$labels) == nrow(d), "entirely unknown support sequence retained in clustering")
stability <- bootstrap_stability(fit, d$patient_id[match(fit$ids, d$hospitalization_id)], 3, 1)
check(nrow(stability) > 0 && all(stability$jaccard >= 0 & stability$jaccard <= 1), "patient bootstrap produces valid stability metrics")
check(identical(rectangular_assignment(matrix(c(9, 8, 1, 7), nrow = 1)), 3L), "fewer local clusters can select any prototype")
check(sum(is.na(rectangular_assignment(matrix(c(1, 2, 3), ncol = 1)))) == 2, "excess local clusters remain explicitly unmatched")
assignment <- tibble::tibble(hospitalization_id = fit$ids, local_cluster = fit$labels)
sig <- cluster_signature(panel, assignment, "joint"); sig$site <- "SyntheticA"
sig2 <- sig; sig2$site <- "SyntheticB"; signatures <- dplyr::bind_rows(sig, sig2)
consensus <- fit_consensus(signatures)
check(all(consensus$mapping$matched), "outcome-free consensus matches complete signatures")
check(all(consensus_loso(signatures, consensus)$held_out_agreement == 1), "leave-one-site-out recovers equivalent synthetic profiles")

options(clif.aggregate_policy='legacy')
masked <- protect_aggregate(tibble::tibble(group = "a", n = c(2, 98), probability = c(.02, .98)), "n", "group", 10, c("n", "probability"))
check(all(is.na(masked$n)) && all(masked$cell_status == "suppressed_partition"), "complementary suppression masks whole partition")
options(clif.aggregate_policy="exact")
check(expect_error(write_share(tibble::tibble(patient_id = "p"), tempdir(), "unsafe", cfg)), "patient IDs rejected from share export")

clinical <- land |> dplyr::left_join(cf |> dplyr::select(hospitalization_id, charlson_noncancer, advanced_cancer), by = "hospitalization_id") |>
  dplyr::left_join(sf, by = "hospitalization_id") |>
  dplyr::left_join(support_summaries(panel), by = "hospitalization_id") |>
  dplyr::left_join(code_status_features(d, tables$code_status, 72, cfg), by = "hospitalization_id") |>
  dplyr::left_join(assignment, by = "hospitalization_id") |>
  dplyr::left_join(consensus$mapping |> dplyr::filter(site == "SyntheticA") |> dplyr::select(local_cluster, phenotype), by = "local_cluster")
aj <- aj_grid(clinical)
repeated <- clinical; repeated$patient_id[2] <- repeated$patient_id[1]
repeated$phenotype[2] <- repeated$phenotype[1]
repeated_aj <- aj_grid(repeated)
check(nrow(repeated_aj) > 0 && all(is.finite(repeated_aj$std_error)), "competing-risk curves retain repeated hospitalization episodes with patient-clustered uncertainty")
check(nrow(aj) > 0 && all(aj$cumulative_incidence >= 0 & aj$cumulative_incidence <= 1), "competing-event cumulative incidence stays within probability bounds")
check(!any(vapply(aj, function(x) is.matrix(x) || is.list(x), logical(1))), "cumulative-incidence exports contain scalar columns")
models <- fit_associations(clinical)
check(nrow(models$metadata) == 4 && setequal(models$metadata$endpoint, cox_endpoints()) && all(models$metadata$fitted_n == nrow(clinical)), "separate death hospice composite live-discharge models report actual fitted denominator")
seq <- fit$sequence[match(clinical$hospitalization_id, fit$ids), , drop = FALSE]
validation <- validate_nested(clinical, seq, cfg)
check(validation$audit$status == "ok", paste("fold-isolated nested validation completes:", validation$audit$status))
check(nrow(validation$metrics) == 12 && nrow(validation$differences) == 8, "validated discrimination calibration Brier and paired increments exported")
pe <- pool_effect(tibble::tibble(site = c("A", "B", "C"), estimate = c(.2, .3, .25), std_error = c(.1, .12, .11)))
check(pe$status == "ok" && is.finite(pe$prediction_low), "REML modified Hartung-Knapp reports prediction interval")
check(expect_error(pool_effect(tibble::tibble(site = c("A", "A"), estimate = c(1, 2), std_error = c(.2, .3)))), "duplicate site effects rejected")

cat("All synthetic regression tests passed. No patient data used.\n")
