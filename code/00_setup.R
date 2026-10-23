required <- c("dplyr", "tidyr", "tibble", "readr", "lubridate", "jsonlite",
              "arrow", "fst", "TraMineR", "cluster", "survival", "comorbidity", "clue", "metafor", "digest")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
if (!all(vapply(required, requireNamespace, logical(1), quietly = TRUE))) stop("Required packages are missing")
message("Setup complete")
