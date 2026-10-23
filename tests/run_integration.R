# Synthetic two-site derivation -> consensus -> model -> pooling round trip.
source("tests/run_tests.R")
scratch <- tempfile("clif_revision_integration_", tmpdir = Sys.getenv("TMPDIR", unset = tempdir()))
dir.create(scratch)
cat("Integration outputs:", scratch, "\n")
input <- file.path(scratch, "tables"); dir.create(input)
tables <- fixture(80L)
# Exercise retention through the full pipeline, not just panel construction:
# h1 has no support observations, h2 has one observation, h3 dies in hour zero.
tables$respiratory_support <- tables$respiratory_support |>
  dplyr::filter(hospitalization_id != "h1") |>
  dplyr::group_by(hospitalization_id) |>
  dplyr::filter(hospitalization_id != "h2" | dplyr::row_number() == 1L) |>
  dplyr::ungroup()
j <- which(tables$hospitalization$hospitalization_id == "h3")
t0_early <- tables$adt$in_dttm[match("h3", tables$adt$hospitalization_id)]
tables$hospitalization$discharge_category[j] <- "Expired"
tables$hospitalization$discharge_dttm[j] <- t0_early + 1800
tables$patient$death_dttm[match(tables$hospitalization$patient_id[j], tables$patient$patient_id)] <- t0_early + 1800
for (nm in names(tables)) readr::write_csv(tables[[nm]], file.path(input, paste0("clif_", nm, ".csv")))
configs <- list(); roots <- list()
for (site in c("SyntheticA", "SyntheticB")) {
  site_cfg <- cfg; site_cfg$site_name <- site; site_cfg$repo <- scratch
  site_cfg$tables_path <- input; site_cfg$analysis$bootstrap_reps <- 0L
  site_cfg$analysis$validation_bootstrap_reps <- 2L
  configs[[site]] <- site_cfg
  prepared <- prepare_revision(site_cfg)
  analyses <- derive_revision(prepared, site_cfg)
  check(length(analyses) == 9L, paste(site, "all default derivations complete"))
  check(all(analyses$h72_arf_onset$data$time_from_landmark > 0), "ARF-anchored outcomes occur after its own landmark")
  check(nrow(analyses$h72_absorbing$data) == 80L, paste(site, "all cohort hospitalizations retained in absorbing analysis"))
  check(nrow(analyses$h72_primary$data) == 79L, paste(site, "all landmark-eligible hospitalizations retained in primary analysis"))
  check(all(c("h1", "h2", "h3") %in% analyses$h72_absorbing$fit$ids), paste(site, "unknown sparse and first-hour death sequences assigned"))
  roots[[site]] <- prepared$root
}
consensus_dir <- file.path(scratch, "consensus")
mapping <- make_consensus_package(file.path(unlist(roots), "share"), consensus_dir)
check(nrow(mapping) > 0 && all(mapping$matched), "two-site consensus package and mapping complete")
check(all(mapping$analysis_role[mapping$analysis_id == "h72_primary"] == "primary") &&
  all(mapping$analysis_role[mapping$analysis_id != "h72_primary"] == "sensitivity"), "primary and sensitivity roles survive central harmonization")
model_dirs <- vapply(names(configs), function(site) {
  model_revision(configs[[site]], file.path(consensus_dir, "consensus_mapping.csv"))
}, character(1))
pooled <- pool_model_packages(model_dirs, file.path(scratch, "pooled"))
check(nrow(pooled) > 0 && any(pooled$status == "ok"), "two-site phenotype effects pooled from fresh model packages")
for (root in roots) {
  exports <- list.files(file.path(root, "share"), recursive = TRUE, full.names = TRUE)
  check(!any(endsWith(exports, ".rds")), "no patient-level RDS in shared package")
  for (f in exports[endsWith(exports, ".csv")]) {
    names <- names(readr::read_csv(f, n_max = 0, show_col_types = FALSE))
    check(!any(c("patient_id", "hospitalization_id", "admission_dttm", "discharge_dttm") %in% names),
          paste("aggregate-only columns", basename(f)))
  }
}
legacy <- file.path(scratch, "legacy.csv")
readr::write_csv(tibble::tibble(n = 1), legacy)
check(expect_error(read_return_csv(legacy)), "legacy return formats rejected")
check(expect_error(make_consensus_package(rep(file.path(roots[[1]], "share"), 2), file.path(scratch, "duplicate"))), "duplicate site derivation packages rejected")
check(expect_error(model_revision(configs[[1]], file.path(consensus_dir, "consensus_mapping.csv"))), "repeat model export cannot overwrite existing consensus package")
cat("Synthetic integration passed. Temporary outputs:", scratch, "\n")
