source("code/lib/load.R")
source("code/lib/central.R")
root <- tempfile("clif_descriptives_"); dir.create(root)
cfg <- list(site_name = "A", analysis = run_defaults())
dirs <- file.path(root, c("A", "B")); for (p in dirs) dir.create(file.path(p, "share"), recursive = TRUE)
mapping <- tibble::tibble(site = c("A", "B"), run_id = cfg$analysis$run_id,
  analysis_id = "h72_primary", local_cluster = 1L, phenotype = "P1")
for (i in 1:2) {
  cfg$site_name <- c("A", "B")[i]
  x <- tibble::tibble(analysis_id = "h72_primary", local_cluster = 1L, variable = "age_years",
    n = 20L, missing_n = 0L, sum = c(1000, 1200)[i], sum_squares = c(51000, 73000)[i], mean = c(50, 60)[i], cell_status = "observed")
  write_share(x, dirs[i], "h72_primary_continuous", cfg)
  cat <- tibble::tibble(analysis_id = "h72_primary", local_cluster = 1L, variable = "death_in_hosp",
    value = c("0", "1"), n = c(15, 5), denominator = 20L, percent = c(75, 25), cell_status = "observed")
  write_share(cat, dirs[i], "h72_primary_categorical", cfg)
}
out <- file.path(root, "pooled"); dir.create(out)
pool_descriptives(file.path(dirs, "share"), mapping, out)
x <- readr::read_csv(file.path(out, "pooled_phenotype_continuous.csv"), show_col_types = FALSE)
y <- readr::read_csv(file.path(out, "pooled_phenotype_categorical.csv"), show_col_types = FALSE)
stopifnot(x$n == 40, x$mean == 55, is.finite(x$sd), all(y$denominator == 40),
          y$n[y$value == 1] == 10, y$percent[y$value == 1] == 25)
cat("PASS: pooled sums, variance, category denominators, and percentages reconcile\n")
cat$cell_status <- "suppressed_partition"; cat$n <- NA_real_; cat$denominator <- NA_real_
cfg$site_name <- "B"; write_share(cat, dirs[2], "h72_primary_categorical", cfg)
pool_descriptives(file.path(dirs, "share"), mapping, out)
y <- readr::read_csv(file.path(out, "pooled_phenotype_categorical.csv"), show_col_types = FALSE)
stopifnot(all(is.na(y$n)), all(is.na(y$denominator)), all(y$status != "complete"))
cat("PASS: suppressed site cells never reconstructed or converted to zero centrally\n")
