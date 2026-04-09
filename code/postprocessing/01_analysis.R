suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

ROOT_DIR <- "."
SITE_DIR <- file.path(ROOT_DIR, "sites")
EXCLUDED_SITES <- c("Hopkins")
RUN_DIRS_ALL <- sort(Sys.glob(file.path(SITE_DIR, "run_*")))
RUN_DATE <- format(Sys.Date(), "%Y%m%d")
OUT_DIR <- file.path(ROOT_DIR, "output", paste0("pooled_", RUN_DATE))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

write_csv_safe <- function(df, path) {
  readr::write_csv(df, path, na = "")
}

save_plot_safe <- function(plot_obj, path, width = 10, height = 7, dpi = 320) {
  ggplot2::ggsave(filename = path, plot = plot_obj, width = width, height = height, dpi = dpi, bg = "white")
}

theme_manuscript <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}

label_resp_support <- function(x) {
  dplyr::recode(
    x,
    "ROOM_AIR" = "Room air",
    "ROOM AIR" = "Room air",
    "LOW_O2" = "Low-flow oxygen",
    "NIV" = "Noninvasive ventilation",
    "IMV" = "Invasive ventilation",
    "OTHER" = "Other support",
    .default = x
  )
}

label_arf <- function(x) {
  dplyr::recode(
    x,
    "NO_ARF" = "No ARF",
    "HYPOX" = "Hypoxemic ARF",
    "HYPER" = "Hypercapnic ARF",
    "MIXED" = "Mixed ARF",
    "NO_ARF_PREDOM" = "Predominantly no ARF",
    "NO_ARF_WITH_HYPOX" = "No ARF with hypoxemia",
    "HYPOX_MIXED_ARF" = "Hypoxemic or mixed ARF",
    "ARF_HYPOX" = "Hypoxemic ARF",
    "ARF_HYPER" = "Hypercapnic ARF",
    "ARF_MIXED" = "Mixed ARF",
    .default = x
  )
}

make_cluster_title <- function(cluster_num, rs_code, arf_code) {
  paste0("Cluster ", cluster_num, ": ", label_resp_support(rs_code), " / ", label_arf(arf_code))
}

label_resp_support_short <- function(x) {
  dplyr::recode(
    x,
    "Room air" = "Room air",
    "Low-flow oxygen" = "Low-flow O2",
    "Noninvasive ventilation" = "NIV",
    "Invasive ventilation" = "IMV",
    "Other support" = "Other",
    .default = x
  )
}

label_arf_short <- function(x) {
  dplyr::recode(
    x,
    "No ARF" = "No ARF",
    "Hypoxemic ARF" = "Hypoxemic ARF",
    "Hypercapnic ARF" = "Hypercapnic ARF",
    "Mixed ARF" = "Mixed ARF",
    "Predominantly no ARF" = "Predom. no ARF",
    "No ARF with hypoxemia" = "No ARF + hypoxemia",
    "Hypoxemic or mixed ARF" = "Hypoxemic/mixed ARF",
    .default = x
  )
}

make_phenotype_label_short <- function(rs_label, arf_label) {
  paste0(label_resp_support_short(rs_label), " / ", label_arf_short(arf_label))
}

pretty_exposure <- function(x) {
  dplyr::recode(
    x,
    "NO2" = "NO₂",
    "PM2.5" = "PM₂.₅",
    "pm25_5y_z" = "PM₂.₅",
    "no2_5y_z" = "NO₂",
    "pm25_mean" = "PM₂.₅",
    "no2_mean" = "NO₂",
    .default = x
  )
}

read_stub_csv <- function(run_dir, stub) {
  hits <- Sys.glob(file.path(run_dir, paste0(stub, "*.csv")))
  if (!length(hits)) return(NULL)
  out <- readr::read_csv(hits[[1]], show_col_types = FALSE, col_types = cols(.default = col_guess()))
  numeric_candidates <- c(
    "traj_cluster_ra", "severity_rank", "n_cluster", "n", "denom", "prop",
    "n_nonmissing", "mean", "sd", "median", "p25", "p75", "estimate", "std_error",
    "statistic", "p_value", "effect", "conf_low", "conf_high", "low", "high",
    "risk_diff", "risk_ratio", "events"
  )
  present_numeric <- intersect(numeric_candidates, names(out))
  if (length(present_numeric)) {
    out <- out %>%
      mutate(across(all_of(present_numeric), ~ suppressWarnings(readr::parse_number(as.character(.x)))))
  }
  out
}

infer_site_run <- function(run_dir) {
  basename(run_dir)
}

infer_site_name <- function(run_dir) {
  nm <- sub("^run_", "", basename(run_dir))
  sub("_[0-9]{8}$", "", nm)
}

RUN_DIRS <- RUN_DIRS_ALL[!infer_site_name(RUN_DIRS_ALL) %in% EXCLUDED_SITES]
stopifnot(length(RUN_DIRS) > 0)

`%||%` <- function(x, y) if (is.null(x)) y else x

harmonize_level <- function(x, variable = NULL) {
  x_chr <- trimws(as.character(x))
  x_low <- tolower(x_chr)
  out <- case_when(
    x_low %in% c("female") ~ "Female",
    x_low %in% c("male") ~ "Male",
    x_low %in% c("unknown", "unk") ~ "Unknown",
    x_low %in% c("hispanic") ~ "Hispanic",
    x_low %in% c("non-hispanic", "non hispanic") ~ "Non-Hispanic",
    x_low %in% c("american indian or alaska native") ~ "American Indian or Alaska Native",
    x_low %in% c("asian") ~ "Asian",
    x_low %in% c("black or african american") ~ "Black or African American",
    x_low %in% c("native hawaiian or other pacific islander") ~ "Native Hawaiian or Other Pacific Islander",
    x_low %in% c("other") ~ "Other",
    x_low %in% c("white") ~ "White",
    x_low %in% c("na", "") ~ "NA",
    TRUE ~ x_chr
  )
  ifelse(is.na(x), NA_character_, out)
}

permute_int <- function(x) {
  if (length(x) <= 1) return(list(x))
  out <- vector("list", factorial(length(x)))
  idx <- 1L
  for (i in seq_along(x)) {
    rest <- x[-i]
    for (p in permute_int(rest)) {
      out[[idx]] <- c(x[[i]], p)
      idx <- idx + 1L
    }
  }
  out
}

best_assignment <- function(cost_mat) {
  k <- nrow(cost_mat)
  perms <- permute_int(seq_len(k))
  best_cost <- Inf
  best_perm <- NULL
  for (perm in perms) {
    cost <- sum(cost_mat[cbind(seq_len(k), perm)])
    if (cost < best_cost) {
      best_cost <- cost
      best_perm <- perm
    }
  }
  tibble(local_idx = seq_len(k), prototype_idx = best_perm, distance = cost_mat[cbind(seq_len(k), best_perm)])
}

pooled_sd <- function(n, mu, s) {
  ok <- !is.na(n) & !is.na(mu) & !is.na(s) & n > 0
  n <- n[ok]
  mu <- mu[ok]
  s <- s[ok]
  if (!length(n)) return(NA_real_)
  if (sum(n) <= 1) return(NA_real_)
  mu_bar <- sum(n * mu) / sum(n)
  ss_within <- sum((n - 1) * s^2)
  ss_between <- sum(n * (mu - mu_bar)^2)
  sqrt((ss_within + ss_between) / (sum(n) - 1))
}

combine_continuous <- function(df, strata_col) {
  groups <- c(strata_col, "variable", "level", "summary_type")
  df %>%
    group_by(across(all_of(groups))) %>%
    group_modify(~ {
      tibble(
        n_nonmissing = sum(.x$n_nonmissing, na.rm = TRUE),
        mean = weighted.mean(.x$mean, .x$n_nonmissing, na.rm = TRUE),
        sd = pooled_sd(.x$n_nonmissing, .x$mean, .x$sd),
        median = NA_real_,
        p25 = NA_real_,
        p75 = NA_real_
      )
    }) %>%
    ungroup()
}

combine_categorical <- function(df, strata_col) {
  groups <- c(strata_col, "variable", "level", "summary_type")
  df %>%
    group_by(across(all_of(groups))) %>%
    summarise(
      n = sum(n, na.rm = TRUE),
      denom = sum(denom, na.rm = TRUE),
      prop = ifelse(sum(denom, na.rm = TRUE) > 0, sum(n, na.rm = TRUE) / sum(denom, na.rm = TRUE), NA_real_),
      .groups = "drop"
    )
}

meta_random_effects <- function(df) {
  df <- df %>%
    filter(!is.na(estimate), !is.na(std_error), std_error > 0)
  k <- nrow(df)
  if (!k) {
    return(tibble(
      k = 0,
      estimate_fixed = NA_real_,
      std_error_fixed = NA_real_,
      conf_low_fixed = NA_real_,
      conf_high_fixed = NA_real_,
      estimate_random = NA_real_,
      std_error_random = NA_real_,
      conf_low_random = NA_real_,
      conf_high_random = NA_real_,
      tau2 = NA_real_,
      q = NA_real_,
      i2 = NA_real_
    ))
  }
  vi <- df$std_error^2
  wi <- 1 / vi
  mu_fixed <- sum(wi * df$estimate) / sum(wi)
  se_fixed <- sqrt(1 / sum(wi))
  q <- sum(wi * (df$estimate - mu_fixed)^2)
  df_q <- k - 1
  c_term <- sum(wi) - (sum(wi^2) / sum(wi))
  tau2 <- if (df_q > 0) max(0, (q - df_q) / c_term) else 0
  wi_re <- 1 / (vi + tau2)
  mu_random <- sum(wi_re * df$estimate) / sum(wi_re)
  se_random <- sqrt(1 / sum(wi_re))
  i2 <- if (df_q > 0 && q > 0) max(0, (q - df_q) / q) else 0
  tibble(
    k = k,
    estimate_fixed = mu_fixed,
    std_error_fixed = se_fixed,
    conf_low_fixed = mu_fixed - 1.96 * se_fixed,
    conf_high_fixed = mu_fixed + 1.96 * se_fixed,
    estimate_random = mu_random,
    std_error_random = se_random,
    conf_low_random = mu_random - 1.96 * se_random,
    conf_high_random = mu_random + 1.96 * se_random,
    tau2 = tau2,
    q = q,
    i2 = i2
  )
}

parse_local_cluster_from_term <- function(term) {
  readr::parse_number(stringr::str_extract(as.character(term), "traj_cluster_ra[0-9]+"))
}

feature_cols <- c(
  "pct_room_air", "pct_low_o2", "pct_niv", "pct_imv", "pct_other",
  "pct_no_arf", "pct_hypox", "pct_hyper", "pct_mixed",
  "pct_vaso", "transition_count"
)

cluster_files <- list(
  centroids = "cluster_centroids_federated",
  severity = "cluster_severity_rank",
  static = "cluster_static_summary",
  rs_signature = "cluster_signature_rs_hourly",
  arf_signature = "cluster_signature_arf_hourly",
  cluster_exposure_mortality = "cluster_exposure_mortality_summary",
  exposure_sofa_corr = "correlation_exposure_sofa",
  overall_cont = "table1_overall_continuous",
  overall_cat = "table1_overall_categorical",
  bycluster_cont = "table1_by_cluster_continuous",
  bycluster_cat = "table1_by_cluster_categorical",
  model_outputs = "model_outputs_standardized",
  model_meta = "model_metadata",
  contrast_no2 = "contrast_no2_low_vs_high_by_cluster",
  contrast_pm25 = "contrast_pm25_low_vs_high_by_cluster",
  run_params = "meta_run_parameters",
  spatial_summary = "spatial_county_cluster_summary",
  spatial_distribution = "spatial_county_cluster_distribution"
)

site_payload <- map(RUN_DIRS, function(run_dir) {
  site_run <- infer_site_run(run_dir)
  site_name <- infer_site_name(run_dir)
  tables <- map(cluster_files, ~ read_stub_csv(run_dir, .x))
  tables$site_run <- site_run
  tables$site_name <- site_name
  tables
})
names(site_payload) <- map_chr(site_payload, "site_run")

centroids <- map_dfr(site_payload, function(x) {
  df <- x$centroids
  if (is.null(df)) return(tibble())
  df %>%
    mutate(
      site_run = x$site_run,
      site_name = x$site_name,
      local_cluster = as.integer(traj_cluster_ra)
    )
})

severity <- map_dfr(site_payload, function(x) {
  df <- x$severity
  if (is.null(df)) return(tibble())
  df %>%
    mutate(
      site_run = x$site_run,
      site_name = x$site_name,
      local_cluster = as.integer(traj_cluster_ra)
    )
})

static_summary <- map_dfr(site_payload, function(x) {
  df <- x$static
  if (is.null(df)) return(tibble())
  df %>%
    mutate(
      site_run = x$site_run,
      site_name = x$site_name,
      local_cluster = as.integer(traj_cluster_ra)
    )
})

stopifnot(nrow(centroids) > 0, nrow(severity) > 0, nrow(static_summary) > 0)

feature_frame <- centroids %>%
  select(site_run, site_name, local_cluster, n_cluster, all_of(feature_cols)) %>%
  mutate(transition_count = as.numeric(scale(transition_count)[, 1]))

feature_matrix <- as.matrix(feature_frame[, feature_cols])
set.seed(20260408)
init_fit <- kmeans(feature_matrix, centers = 5, nstart = 100)
prototypes <- init_fit$centers

assignments <- NULL
for (iter in seq_len(25)) {
  assignments_new <- feature_frame %>%
    group_by(site_run, site_name) %>%
    group_modify(~ {
      cost <- as.matrix(dist(rbind(as.matrix(.x[, feature_cols]), prototypes)))[seq_len(nrow(.x)), nrow(.x) + seq_len(nrow(prototypes))]
      best_assignment(cost) %>%
        mutate(local_cluster = .x$local_cluster[local_idx])
    }) %>%
    ungroup() %>%
    select(site_run, site_name, local_cluster, prototype_idx, distance)

  prototype_updates <- feature_frame %>%
    left_join(assignments_new, by = c("site_run", "site_name", "local_cluster")) %>%
    group_by(prototype_idx) %>%
    summarise(
      across(all_of(feature_cols), ~ weighted.mean(.x, n_cluster, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(prototype_idx)

  prototypes_new <- as.matrix(prototype_updates[, feature_cols])
  rownames(prototypes_new) <- prototype_updates$prototype_idx

  stable <- !is.null(assignments) &&
    identical(assignments_new %>% arrange(site_run, local_cluster) %>% pull(prototype_idx),
              assignments %>% arrange(site_run, local_cluster) %>% pull(prototype_idx))

  assignments <- assignments_new
  prototypes <- prototypes_new
  if (stable) break
}

cluster_map_raw <- assignments %>%
  mutate(prototype_idx = as.integer(prototype_idx))

pooled_cluster_summary <- centroids %>%
  left_join(cluster_map_raw, by = c("site_run", "site_name", "local_cluster")) %>%
  left_join(
    static_summary %>%
      select(site_run, local_cluster, n, death_pct, icu_los_median, imv72_pct, sofa_total_mean),
    by = c("site_run", "local_cluster")
  ) %>%
  group_by(prototype_idx) %>%
  summarise(
    n_sites = n_distinct(site_run),
    n_patients = sum(n_cluster, na.rm = TRUE),
    across(all_of(feature_cols), ~ weighted.mean(.x, n_cluster, na.rm = TRUE)),
    death_pct = weighted.mean(death_pct, n_cluster, na.rm = TRUE),
    icu_los_median_proxy = weighted.mean(icu_los_median, n_cluster, na.rm = TRUE),
    imv72_pct = weighted.mean(imv72_pct, n_cluster, na.rm = TRUE),
    sofa_total_mean = weighted.mean(sofa_total_mean, n_cluster, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    respiratory_burden = scale(pct_imv + pct_niv + pct_vaso + pct_hypox + pct_mixed)[, 1],
    clinical_burden = scale(sofa_total_mean + death_pct + imv72_pct + icu_los_median_proxy / 100)[, 1],
    severity_score = respiratory_burden + clinical_burden,
    dominant_rs = c("ROOM_AIR", "LOW_O2", "NIV", "IMV", "OTHER")[max.col(select(., pct_room_air, pct_low_o2, pct_niv, pct_imv, pct_other), ties.method = "first")],
    dominant_arf = c("NO_ARF", "HYPOX", "HYPER", "MIXED")[max.col(select(., pct_no_arf, pct_hypox, pct_hyper, pct_mixed), ties.method = "first")]
  ) %>%
  arrange(severity_score) %>%
  mutate(
    pooled_cluster = row_number(),
    pooled_severity_rank = row_number(),
    is_worst_cluster = pooled_cluster == max(pooled_cluster),
    arf_label = case_when(
      dominant_rs == "IMV" & (pct_hypox + pct_mixed) >= 0.45 ~ "HYPOX_MIXED_ARF",
      dominant_rs == "IMV" & pct_no_arf >= 0.75 ~ "NO_ARF_PREDOM",
      dominant_arf == "NO_ARF" & pct_hypox >= 0.2 ~ "NO_ARF_WITH_HYPOX",
      TRUE ~ dominant_arf
    ),
    pooled_cluster_label = make_cluster_title(pooled_cluster, dominant_rs, arf_label),
    phenotype_label = paste0(label_resp_support(dominant_rs), " / ", label_arf(arf_label)),
    phenotype_label_short = make_phenotype_label_short(label_resp_support(dominant_rs), label_arf(arf_label))
  )

prototype_key <- pooled_cluster_summary %>%
  select(prototype_idx, pooled_cluster, pooled_severity_rank, arf_label, dominant_rs, pooled_cluster_label, phenotype_label, phenotype_label_short, is_worst_cluster, severity_score)

cluster_map <- cluster_map_raw %>%
  left_join(prototype_key, by = "prototype_idx") %>%
  left_join(
    severity %>%
      select(site_run, local_cluster, local_severity_rank = severity_rank),
    by = c("site_run", "local_cluster")
  ) %>%
  left_join(
    static_summary %>%
      select(site_run, local_cluster, n, death_pct, imv72_pct, sofa_total_mean, icu_los_median),
    by = c("site_run", "local_cluster")
  ) %>%
  arrange(site_run, pooled_cluster)

write_csv_safe(cluster_map, file.path(OUT_DIR, "pooled_cluster_mapping.csv"))
write_csv_safe(pooled_cluster_summary, file.path(OUT_DIR, "pooled_cluster_prototypes.csv"))

overall_cont <- map_dfr(site_payload, function(x) {
  df <- x$overall_cont
  if (is.null(df)) return(tibble())
  df %>% mutate(site_run = x$site_run, site_name = x$site_name)
}) %>%
  mutate(stratum = "Overall")

overall_cat <- map_dfr(site_payload, function(x) {
  df <- x$overall_cat
  if (is.null(df)) return(tibble())
  df %>%
    mutate(
      level = harmonize_level(level, variable),
      site_run = x$site_run,
      site_name = x$site_name
    )
}) %>%
  mutate(stratum = "Overall")

pooled_overall_cont <- overall_cont %>%
  combine_continuous("stratum") %>%
  mutate(site_name = "pooled") %>%
  relocate(site_name, stratum, variable, level, summary_type)

pooled_overall_cat_base <- overall_cat %>%
  filter(variable != "traj_cluster_ra") %>%
  combine_categorical("stratum") %>%
  mutate(site_name = "pooled")

pooled_cluster_dist <- cluster_map %>%
  transmute(
    stratum = "Overall",
    variable = "pooled_cluster",
    level = as.character(pooled_cluster),
    summary_type = "categorical",
    n = n,
    denom = sum(n, na.rm = TRUE),
    prop = n / denom
  ) %>%
  select(stratum, variable, level, summary_type, n, denom, prop) %>%
  mutate(site_name = "pooled")

pooled_overall_cat <- bind_rows(pooled_overall_cat_base, pooled_cluster_dist) %>%
  arrange(variable, level) %>%
  relocate(site_name, stratum, variable, level, summary_type)

bycluster_cont <- map_dfr(site_payload, function(x) {
  df <- x$bycluster_cont
  if (is.null(df)) return(tibble())
  df %>% mutate(site_run = x$site_run, site_name = x$site_name)
})

bycluster_cat <- map_dfr(site_payload, function(x) {
  df <- x$bycluster_cat
  if (is.null(df)) return(tibble())
  df %>%
    mutate(
      level = harmonize_level(level, variable),
      site_run = x$site_run,
      site_name = x$site_name
    )
})

pooled_bycluster_cont <- bycluster_cont %>%
  mutate(local_cluster = suppressWarnings(as.integer(stratum))) %>%
  filter(!is.na(local_cluster)) %>%
  left_join(cluster_map %>% select(site_run, local_cluster, pooled_cluster, pooled_cluster_label), by = c("site_run", "local_cluster")) %>%
  filter(!is.na(pooled_cluster)) %>%
  mutate(stratum = as.character(pooled_cluster)) %>%
  combine_continuous("stratum") %>%
  mutate(site_name = "pooled") %>%
  left_join(distinct(cluster_map, pooled_cluster, pooled_cluster_label) %>% mutate(stratum = as.character(pooled_cluster)), by = "stratum") %>%
  relocate(site_name, stratum, pooled_cluster_label, variable, level, summary_type)

pooled_bycluster_cat <- bycluster_cat %>%
  mutate(local_cluster = suppressWarnings(as.integer(stratum))) %>%
  filter(!is.na(local_cluster)) %>%
  left_join(cluster_map %>% select(site_run, local_cluster, pooled_cluster, pooled_cluster_label), by = c("site_run", "local_cluster")) %>%
  filter(!is.na(pooled_cluster)) %>%
  mutate(stratum = as.character(pooled_cluster)) %>%
  combine_categorical("stratum") %>%
  mutate(site_name = "pooled") %>%
  left_join(distinct(cluster_map, pooled_cluster, pooled_cluster_label) %>% mutate(stratum = as.character(pooled_cluster)), by = "stratum") %>%
  relocate(site_name, stratum, pooled_cluster_label, variable, level, summary_type)

write_csv_safe(pooled_overall_cont, file.path(OUT_DIR, "pooled_table1_overall_continuous.csv"))
write_csv_safe(pooled_overall_cat, file.path(OUT_DIR, "pooled_table1_overall_categorical.csv"))
write_csv_safe(pooled_bycluster_cont, file.path(OUT_DIR, "pooled_table1_by_cluster_continuous.csv"))
write_csv_safe(pooled_bycluster_cat, file.path(OUT_DIR, "pooled_table1_by_cluster_categorical.csv"))

fmt_n_pct <- function(n, denom) {
  ifelse(is.na(n) | is.na(denom) | denom == 0, NA_character_, sprintf("%s (%.1f%%)", scales::comma(n), 100 * n / denom))
}

fmt_mean_sd <- function(mean, sd) {
  ifelse(is.na(mean), NA_character_, sprintf("%.1f (%.1f)", mean, sd))
}

fmt_effect <- function(est, lo, hi, digits = 2) {
  ifelse(is.na(est), NA_character_, sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"), est, lo, hi))
}

cluster_ns <- cluster_map %>%
  group_by(pooled_cluster, pooled_cluster_label) %>%
  summarise(n = sum(n, na.rm = TRUE), .groups = "drop")

overall_n <- sum(cluster_ns$n, na.rm = TRUE)

cont_labels <- c(
  age_years = "Age, mean (SD)",
  charlson_score = "Charlson score, mean (SD)",
  sofa_total = "SOFA total, mean (SD)",
  pm25_5y = "PM₂.₅ 5-year exposure, mean (SD)",
  no2_5y = "NO₂ 5-year exposure, mean (SD)",
  icu_los_hours = "ICU LOS (hours), mean (SD)",
  imv_hours_72h = "IMV hours in first 72h, mean (SD)"
)

cat_specs <- tribble(
  ~variable, ~level, ~label,
  "sex_category", "Male", "Male sex, n (%)",
  "race_category", "White", "White race, n (%)",
  "race_category", "Black or African American", "Black race, n (%)",
  "ethnicity_category", "Hispanic", "Hispanic ethnicity, n (%)",
  "advanced_cancer_any_poa", "1", "Advanced cancer, n (%)",
  "death_in_hosp", "1", "In-hospital death, n (%)",
  "hospice_discharge", "1", "Hospice discharge, n (%)",
  "death_or_hospice", "1", "Death or hospice, n (%)",
  "imv_72h_any", "1", "Any IMV in first 72h, n (%)"
)

clean_overall_cont <- pooled_overall_cont %>%
  filter(variable %in% names(cont_labels), stratum == "Overall") %>%
  transmute(label = unname(cont_labels[variable]), Overall = fmt_mean_sd(mean, sd))

clean_cluster_cont <- pooled_bycluster_cont %>%
  filter(variable %in% names(cont_labels)) %>%
  transmute(
    pooled_cluster,
    label = unname(cont_labels[variable]),
    value = fmt_mean_sd(mean, sd)
  )

clean_overall_cat <- pooled_overall_cat %>%
  inner_join(cat_specs, by = c("variable", "level")) %>%
  transmute(label, Overall = fmt_n_pct(n, denom))

clean_cluster_cat <- pooled_bycluster_cat %>%
  inner_join(cat_specs, by = c("variable", "level")) %>%
  transmute(
    pooled_cluster,
    label,
    value = fmt_n_pct(n, denom)
  )

cluster_n_row <- cluster_ns %>%
  transmute(pooled_cluster, label = "N", value = scales::comma(n))

overall_n_row <- tibble(label = "N", Overall = scales::comma(overall_n))

manuscript_table1 <- bind_rows(
  overall_n_row,
  clean_overall_cont,
  clean_overall_cat
) %>%
  full_join(
    bind_rows(cluster_n_row, clean_cluster_cont, clean_cluster_cat) %>%
      mutate(cluster_col = paste0("Cluster ", pooled_cluster)) %>%
      select(-pooled_cluster) %>%
      pivot_wider(names_from = cluster_col, values_from = value),
    by = "label"
  ) %>%
  mutate(
    sort_order = match(label, c(
      "N",
      unname(cont_labels),
      cat_specs$label
    ))
  ) %>%
  arrange(sort_order) %>%
  select(-sort_order)

write_csv_safe(manuscript_table1, file.path(OUT_DIR, "manuscript_table1.csv"))

model_outputs <- map_dfr(site_payload, function(x) {
  df <- x$model_outputs
  if (is.null(df)) return(tibble())
  df %>% mutate(site_run = x$site_run, site_name = x$site_name)
})

stable_terms <- tribble(
  ~model, ~term_group,
  "polr_severity_rank", "trajectory_severity",
  "glm_mortality", "mortality",
  "glm_los_quasipoisson", "icu_los",
  "cox72_landmark", "post72_survival",
  "glm_total_no2", "mortality_total_no2",
  "glm_direct_no2_severity", "mortality_direct_no2_plus_severity",
  "glm_direct_no2_cluster", "mortality_direct_no2_plus_cluster"
) %>%
  crossing(term = c("pm25_5y_z", "no2_5y_z")) %>%
  filter(!(str_detect(model, "no2") & term == "pm25_5y_z"))

pooled_meta <- model_outputs %>%
  inner_join(stable_terms, by = c("model", "term")) %>%
  group_by(term_group, model, term, effect_type) %>%
  group_modify(~ {
    effect_type_val <- .y$effect_type[[1]] %||% NA_character_
    meta_random_effects(.x) %>%
      mutate(
        pooled_effect_fixed = ifelse(!is.na(effect_type_val) && effect_type_val != "", exp(estimate_fixed), NA_real_),
        pooled_effect_random = ifelse(!is.na(effect_type_val) && effect_type_val != "", exp(estimate_random), NA_real_),
        pooled_conf_low_random = ifelse(!is.na(effect_type_val) && effect_type_val != "", exp(conf_low_random), NA_real_),
        pooled_conf_high_random = ifelse(!is.na(effect_type_val) && effect_type_val != "", exp(conf_high_random), NA_real_),
        site_included = paste(sort(unique(.x$site_name)), collapse = "; ")
      )
  }) %>%
  ungroup()

write_csv_safe(pooled_meta, file.path(OUT_DIR, "pooled_meta_analysis.csv"))

contrast_no2 <- map_dfr(site_payload, function(x) {
  df <- x$contrast_no2
  if (is.null(df)) return(tibble())
  df %>%
    mutate(
      site_run = x$site_run,
      site_name = x$site_name,
      local_cluster = as.integer(traj_cluster_ra)
    )
})

contrast_pm25 <- map_dfr(site_payload, function(x) {
  df <- x$contrast_pm25
  if (is.null(df)) return(tibble())
  df %>%
    mutate(
      site_run = x$site_run,
      site_name = x$site_name,
      local_cluster = as.integer(traj_cluster_ra)
    )
})

summarise_modeled_contrasts <- function(df, exposure_name) {
  df %>%
    left_join(cluster_map %>% select(site_run, local_cluster, pooled_cluster, pooled_cluster_label, n), by = c("site_run", "local_cluster")) %>%
    filter(!is.na(pooled_cluster)) %>%
    group_by(pooled_cluster, pooled_cluster_label) %>%
    summarise(
      exposure = exposure_name,
      n_sites = n_distinct(site_run),
      n_patients = sum(n, na.rm = TRUE),
      risk_diff_weighted = weighted.mean(risk_diff, n, na.rm = TRUE),
      risk_ratio_weighted = exp(weighted.mean(log(risk_ratio), n, na.rm = TRUE)),
      low_prob_weighted = weighted.mean(low, n, na.rm = TRUE),
      high_prob_weighted = weighted.mean(high, n, na.rm = TRUE),
      min_site_risk_diff = min(risk_diff, na.rm = TRUE),
      max_site_risk_diff = max(risk_diff, na.rm = TRUE),
      .groups = "drop"
    )
}

pooled_modeled_contrasts <- bind_rows(
  summarise_modeled_contrasts(contrast_pm25, "pm25_5y_z"),
  summarise_modeled_contrasts(contrast_no2, "no2_5y_z")
) %>%
  arrange(exposure, pooled_cluster)

write_csv_safe(pooled_modeled_contrasts, file.path(OUT_DIR, "pooled_modeled_cluster_contrasts.csv"))

meta_plot_df <- pooled_meta %>%
  mutate(
    exposure = pretty_exposure(recode(term, pm25_5y_z = "PM2.5", no2_5y_z = "NO2")),
    outcome = recode(
      term_group,
      mortality = "Mortality",
      icu_los = "ICU LOS",
      post72_survival = "Post-72h survival",
      trajectory_severity = "Trajectory severity",
      mortality_total_no2 = "NO2 total effect",
      mortality_direct_no2_plus_severity = "NO2 direct effect + severity",
      mortality_direct_no2_plus_cluster = "NO2 direct effect + cluster"
    ),
    effect_display = fmt_effect(pooled_effect_random, pooled_conf_low_random, pooled_conf_high_random)
  ) %>%
  filter(!is.na(pooled_effect_random))

p_meta <- ggplot(
  meta_plot_df,
  aes(x = pooled_effect_random, y = forcats::fct_rev(interaction(outcome, exposure, sep = " | ")), xmin = pooled_conf_low_random, xmax = pooled_conf_high_random, color = exposure)
) +
  geom_vline(xintercept = 1, linetype = 2, color = "gray45") +
  geom_errorbar(width = 0.18, linewidth = 0.8, orientation = "y") +
  geom_point(size = 2.6) +
  scale_x_log10() +
  labs(
    title = "Pooled Exposure Associations Across Outcomes",
    x = "Pooled effect estimate (random-effects, log scale)",
    y = NULL,
    color = NULL
  ) +
  theme_manuscript()

save_plot_safe(p_meta, file.path(OUT_DIR, "forest_overall_meta_effects.png"), width = 11, height = 7)

interaction_terms <- model_outputs %>%
  filter(model %in% c("glm_interaction_no2", "glm_interaction_pm25")) %>%
  filter(str_detect(term, ":traj_cluster_ra|traj_cluster_ra[0-9]+:")) %>%
  mutate(
    exposure = ifelse(model == "glm_interaction_no2", "NO2", "PM2.5"),
    local_cluster = parse_local_cluster_from_term(term)
  ) %>%
  left_join(
    cluster_map %>% select(site_run, local_cluster, pooled_cluster, pooled_cluster_label),
    by = c("site_run", "local_cluster")
  ) %>%
  filter(!is.na(pooled_cluster)) %>%
  mutate(
    recoded_term = paste0(exposure, " x ", pooled_cluster_label),
    site_effect = effect,
    site_conf_low = conf_low,
    site_conf_high = conf_high
  )

pooled_interaction_meta <- interaction_terms %>%
  group_by(model, exposure, pooled_cluster, pooled_cluster_label, recoded_term, effect_type) %>%
  group_modify(~ {
    meta_random_effects(.x) %>%
      mutate(
        pooled_effect_random = exp(estimate_random),
        pooled_conf_low_random = exp(conf_low_random),
        pooled_conf_high_random = exp(conf_high_random),
        site_included = paste(sort(unique(.x$site_name)), collapse = "; ")
      )
  }) %>%
  ungroup() %>%
  arrange(exposure, pooled_cluster)

write_csv_safe(pooled_interaction_meta, file.path(OUT_DIR, "pooled_interaction_meta_analysis.csv"))

interaction_reference_rows <- cluster_map %>%
  filter(local_cluster == 1) %>%
  select(site_run, site_name, pooled_cluster, pooled_cluster_label) %>%
  tidyr::crossing(exposure = c("NO2", "PM2.5")) %>%
  mutate(
    model = ifelse(exposure == "NO2", "glm_interaction_no2", "glm_interaction_pm25"),
    term = "reference",
    estimate = 0,
    std_error = NA_real_,
    effect = 1,
    conf_low = NA_real_,
    conf_high = NA_real_,
    is_reference = TRUE
  )

interaction_plot_terms <- bind_rows(
  interaction_terms %>%
    transmute(
      site_run, site_name, exposure, pooled_cluster, pooled_cluster_label,
      model, term, estimate, std_error, effect = site_effect, conf_low = site_conf_low, conf_high = site_conf_high,
      is_reference = FALSE
    ),
  interaction_reference_rows
)

write_csv_safe(interaction_plot_terms, file.path(OUT_DIR, "interaction_plot_terms.csv"))

interaction_site_order <- c(sort(unique(cluster_map$site_name)), "Pooled random-effects")
interaction_site_order_plot <- rev(interaction_site_order)

plot_interaction_forest <- function(site_df, pooled_df, exposure_label) {
  exposure_display <- pretty_exposure(exposure_label)
  site_plot_df <- site_df %>%
    filter(exposure == exposure_label) %>%
    transmute(
      pooled_cluster_label,
      site_name,
      effect,
      conf_low,
      conf_high,
      n = 1,
      pooled_flag = FALSE,
      is_reference = is_reference
    )

  pooled_plot_df <- pooled_df %>%
    filter(exposure == exposure_label) %>%
    transmute(
      pooled_cluster_label,
      site_name = "Pooled random-effects",
      effect = pooled_effect_random,
      conf_low = pooled_conf_low_random,
      conf_high = pooled_conf_high_random,
      n = k,
      pooled_flag = TRUE,
      is_reference = FALSE
    )

  plot_df <- bind_rows(site_plot_df, pooled_plot_df) %>%
    mutate(
      site_name = factor(site_name, levels = interaction_site_order_plot),
      point_size = ifelse(pooled_flag, 3, 2),
      label_text = case_when(
        is_reference ~ "Ref",
        !is.na(conf_low) & !is.na(conf_high) ~ sprintf("%.2f (%.2f, %.2f)", effect, conf_low, conf_high),
        TRUE ~ sprintf("%.2f", effect)
      ),
      label_x = case_when(
        !is.na(conf_high) ~ conf_high * 1.06,
        TRUE ~ effect * 1.06
      )
    )

  ggplot(plot_df, aes(x = effect, y = site_name, xmin = conf_low, xmax = conf_high, color = pooled_flag)) +
    geom_vline(xintercept = 1, linetype = 2, color = "gray50") +
    geom_errorbar(data = ~ dplyr::filter(.x, !is.na(conf_low), !is.na(conf_high)), width = 0.18, linewidth = 0.7, orientation = "y", show.legend = FALSE) +
    geom_point(aes(size = point_size, shape = is_reference), show.legend = FALSE) +
    geom_text(aes(x = label_x, label = label_text), hjust = 0, size = 2.6, show.legend = FALSE, color = "gray15") +
    scale_x_log10() +
    scale_y_discrete(limits = interaction_site_order_plot, drop = FALSE) +
    scale_color_manual(values = c(`FALSE` = "#6a8ba3", `TRUE` = "#b33a3a")) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1)) +
    scale_size_identity() +
    facet_wrap(~ pooled_cluster_label, ncol = 2, scales = "fixed") +
    labs(
      title = paste0(exposure_display, " Interaction Terms by Consensus Trajectory Cluster"),
      x = "Interaction OR for mortality (log scale)",
      y = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_manuscript() +
    theme(
      plot.margin = margin(10, 90, 10, 10),
      panel.spacing.x = unit(1.2, "lines")
    )
}

save_plot_safe(
  plot_interaction_forest(interaction_plot_terms, pooled_interaction_meta, "NO2"),
  file.path(OUT_DIR, "forest_cluster_interaction_no2.png"),
  width = 12,
  height = 9
)

save_plot_safe(
  plot_interaction_forest(interaction_plot_terms, pooled_interaction_meta, "PM2.5"),
  file.path(OUT_DIR, "forest_cluster_interaction_pm25.png"),
  width = 12,
  height = 9
)

rs_signature <- map_dfr(site_payload, function(x) {
  df <- x$rs_signature
  if (is.null(df)) return(tibble())
  df %>%
    mutate(
      site_run = x$site_run,
      site_name = x$site_name,
      local_cluster = as.integer(traj_cluster_ra),
      rs = label_resp_support(rs)
    )
})

arf_signature <- map_dfr(site_payload, function(x) {
  df <- x$arf_signature
  if (is.null(df)) return(tibble())
  df %>%
    mutate(
      site_run = x$site_run,
      site_name = x$site_name,
      local_cluster = as.integer(traj_cluster_ra),
      arf = label_arf(arf)
    )
})

rs_signature_pooled <- rs_signature %>%
  left_join(cluster_map %>% select(site_run, local_cluster, pooled_cluster, pooled_cluster_label), by = c("site_run", "local_cluster")) %>%
  filter(!is.na(pooled_cluster)) %>%
  group_by(pooled_cluster, pooled_cluster_label, h, rs) %>%
  summarise(n = sum(n, na.rm = TRUE), .groups = "drop_last") %>%
  mutate(prop = n / sum(n, na.rm = TRUE)) %>%
  ungroup()

arf_signature_pooled <- arf_signature %>%
  left_join(cluster_map %>% select(site_run, local_cluster, pooled_cluster, pooled_cluster_label), by = c("site_run", "local_cluster")) %>%
  filter(!is.na(pooled_cluster)) %>%
  group_by(pooled_cluster, pooled_cluster_label, h, arf) %>%
  summarise(n = sum(n, na.rm = TRUE), .groups = "drop_last") %>%
  mutate(prop = n / sum(n, na.rm = TRUE)) %>%
  ungroup()

write_csv_safe(rs_signature_pooled, file.path(OUT_DIR, "pooled_signature_resp_support_hourly.csv"))
write_csv_safe(arf_signature_pooled, file.path(OUT_DIR, "pooled_signature_arf_hourly.csv"))

rs_palette <- c(
  "Room air" = "#d8e4bc",
  "Low-flow oxygen" = "#f4d06f",
  "Noninvasive ventilation" = "#7fb7be",
  "Invasive ventilation" = "#d96c75",
  "Other support" = "#9d9d9d"
)

arf_palette <- c(
  "No ARF" = "#d8e4bc",
  "Hypoxemic ARF" = "#5ca4a9",
  "Hypercapnic ARF" = "#f6aa1c",
  "Mixed ARF" = "#b33a3a"
)

p_rs_sig <- ggplot(rs_signature_pooled, aes(x = h, y = prop, fill = rs)) +
  geom_area(color = "white", linewidth = 0.15, alpha = 0.98) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = rs_palette, name = "Respiratory support") +
  facet_wrap(~ pooled_cluster_label, ncol = 1) +
  labs(
    title = "Pooled Consortium Respiratory Support Trajectories",
    x = "Hours from ICU admission",
    y = "Patients in cluster-state"
  ) +
  theme_manuscript() +
  theme(panel.grid.major.x = element_line(color = "gray90"))

p_arf_sig <- ggplot(arf_signature_pooled, aes(x = h, y = prop, fill = arf)) +
  geom_area(color = "white", linewidth = 0.15, alpha = 0.98) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = arf_palette, name = "ARF subtype") +
  facet_wrap(~ pooled_cluster_label, ncol = 1) +
  labs(
    title = "Pooled Consortium Acute Respiratory Failure Trajectories",
    x = "Hours from ICU admission",
    y = "Patients in cluster-state"
  ) +
  theme_manuscript() +
  theme(panel.grid.major.x = element_line(color = "gray90"))

save_plot_safe(p_rs_sig, file.path(OUT_DIR, "figure_pooled_resp_support_trajectories.png"), width = 11, height = 14)
save_plot_safe(p_arf_sig, file.path(OUT_DIR, "figure_pooled_arf_trajectories.png"), width = 11, height = 14)

severe_arf_df <- arf_signature_pooled %>%
  mutate(severe_arf = arf %in% c("Hypoxemic ARF", "Mixed ARF")) %>%
  group_by(pooled_cluster, pooled_cluster_label, h) %>%
  summarise(prop = sum(prop[severe_arf], na.rm = TRUE), .groups = "drop")

p_severe_arf <- severe_arf_df %>%
  left_join(distinct(prototype_key, pooled_cluster, phenotype_label_short), by = "pooled_cluster") %>%
  ggplot(aes(x = h, y = prop, color = stringr::str_wrap(phenotype_label_short, width = 20))) +
  geom_line(linewidth = 1.1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Severe ARF Burden Over Time by Consensus Trajectory Cluster",
    x = "Hours from ICU admission",
    y = "Share with severe ARF",
    color = "Trajectory phenotype"
  ) +
  theme_manuscript() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(face = "bold")
  )

save_plot_safe(p_severe_arf, file.path(OUT_DIR, "figure_pooled_severe_arf_burden.png"), width = 11, height = 6)

cluster_exposure_mortality <- map_dfr(site_payload, function(x) {
  df <- x$cluster_exposure_mortality
  if (is.null(df)) return(tibble())
  df %>%
    mutate(
      site_run = x$site_run,
      site_name = x$site_name,
      local_cluster = suppressWarnings(as.integer(traj_cluster_ra))
    ) %>%
    filter(!is.na(local_cluster))
})

exposure_sofa_corr <- map_dfr(site_payload, function(x) {
  df <- x$exposure_sofa_corr
  if (is.null(df)) return(tibble())
  df %>%
    mutate(site_run = x$site_run, site_name = x$site_name)
})

cluster_exposure_outcomes <- cluster_exposure_mortality %>%
  left_join(cluster_map %>% select(site_run, local_cluster, pooled_cluster, pooled_cluster_label, n), by = c("site_run", "local_cluster")) %>%
  filter(!is.na(pooled_cluster)) %>%
  group_by(pooled_cluster, pooled_cluster_label) %>%
  summarise(
    n_sites = n_distinct(site_name),
    n_patients = sum(n, na.rm = TRUE),
    mortality = weighted.mean(mortality, n, na.rm = TRUE),
    pm25_mean = weighted.mean(pm25, n, na.rm = TRUE),
    no2_mean = weighted.mean(no2, n, na.rm = TRUE),
    .groups = "drop"
  )

write_csv_safe(cluster_exposure_outcomes, file.path(OUT_DIR, "publication_table_cluster_exposure_outcomes.csv"))

pollution_outcome_table <- bind_rows(
  pooled_meta %>%
    mutate(
      table_section = "Overall pooled associations",
      label = case_when(
        term_group == "mortality" & term == "pm25_5y_z" ~ "PM₂.₅ -> death or hospice",
        term_group == "mortality" & term == "no2_5y_z" ~ "NO₂ -> death or hospice",
        term_group == "icu_los" & term == "pm25_5y_z" ~ "PM₂.₅ -> ICU LOS",
        term_group == "icu_los" & term == "no2_5y_z" ~ "NO₂ -> ICU LOS",
        term_group == "post72_survival" & term == "pm25_5y_z" ~ "PM₂.₅ -> post-72h mortality hazard",
        term_group == "post72_survival" & term == "no2_5y_z" ~ "NO₂ -> post-72h mortality hazard",
        term_group == "trajectory_severity" & term == "pm25_5y_z" ~ "PM₂.₅ -> higher trajectory severity",
        term_group == "trajectory_severity" & term == "no2_5y_z" ~ "NO₂ -> higher trajectory severity",
        TRUE ~ NA_character_
      ),
      estimate_ci = fmt_effect(pooled_effect_random, pooled_conf_low_random, pooled_conf_high_random),
      metric = effect_type
    ) %>%
    filter(!is.na(label)) %>%
    transmute(section = table_section, label, metric, estimate_ci, k, site_included),
  pooled_interaction_meta %>%
    mutate(
      section = "Exposure-by-cluster interaction terms",
      label = gsub("^NO2", "NO₂", gsub("^PM2.5", "PM₂.₅", recoded_term)),
      metric = effect_type,
      estimate_ci = fmt_effect(pooled_effect_random, pooled_conf_low_random, pooled_conf_high_random)
    ) %>%
    transmute(section, label, metric, estimate_ci, k, site_included)
)

write_csv_safe(pollution_outcome_table, file.path(OUT_DIR, "publication_table_air_pollution_results.csv"))

interaction_heatmap_df <- pooled_interaction_meta %>%
  transmute(
    exposure = pretty_exposure(exposure),
    pooled_cluster_label,
    pooled_or = pooled_effect_random,
    conf_low = pooled_conf_low_random,
    conf_high = pooled_conf_high_random
  )

p_interaction_heatmap <- ggplot(interaction_heatmap_df, aes(x = exposure, y = forcats::fct_rev(pooled_cluster_label), fill = pooled_or)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.2f", pooled_or)), size = 3.5) +
  scale_fill_gradient2(
    low = "#2c7fb8",
    mid = "#f7f7f7",
    high = "#b2182b",
    midpoint = 1,
    name = "Pooled\ninteraction OR"
  ) +
  labs(
    title = "Air Pollution Interaction Signals Across Consensus Trajectory Phenotypes",
    x = NULL,
    y = NULL
  ) +
  theme_manuscript()

save_plot_safe(p_interaction_heatmap, file.path(OUT_DIR, "figure_air_pollution_interaction_heatmap.png"), width = 9, height = 5.5)

cluster_exposure_display <- cluster_exposure_outcomes %>%
  left_join(distinct(prototype_key, pooled_cluster, phenotype_label_short), by = "pooled_cluster") %>%
  transmute(
    pooled_cluster,
    phenotype = paste0("C", pooled_cluster, ": ", phenotype_label_short),
    Mortality = mortality * 100,
    `PM₂.₅` = pm25_mean,
    `NO₂` = no2_mean
  ) %>%
  pivot_longer(cols = c(Mortality, `PM₂.₅`, `NO₂`), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = factor(metric, levels = c("Mortality", "PM₂.₅", "NO₂")),
    label = case_when(
      metric == "Mortality" ~ sprintf("%.1f%%", value),
      TRUE ~ sprintf("%.2f", value)
    )
  )

p_cluster_exposure <- ggplot(
  cluster_exposure_display,
  aes(x = metric, y = forcats::fct_rev(phenotype), fill = value)
) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 3.6) +
  facet_wrap(~ metric, scales = "free_x") +
  scale_fill_gradientn(
    colors = c("#f7fbff", "#9ecae1", "#3182bd", "#08519c"),
    name = "Value"
  ) +
  labs(
    title = "Cluster-Level Air Pollution Exposure and Mortality Summary",
    x = NULL,
    y = NULL
  ) +
  theme_manuscript() +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

save_plot_safe(p_cluster_exposure, file.path(OUT_DIR, "figure_cluster_exposure_mortality.png"), width = 10, height = 5.8)

sofa_corr_summary <- exposure_sofa_corr %>%
  mutate(
    variable_x = recode(variable_x, pm25_5y = "PM₂.₅", no2_5y = "NO₂")
  ) %>%
  group_by(variable_x, variable_y) %>%
  summarise(
    mean_correlation = mean(estimate, na.rm = TRUE),
    median_correlation = median(estimate, na.rm = TRUE),
    n_sites = n_distinct(site_name),
    .groups = "drop"
  )

write_csv_safe(sofa_corr_summary, file.path(OUT_DIR, "publication_table_exposure_sofa_correlations.csv"))

spatial_distribution <- map_dfr(site_payload, function(x) {
  df <- x$spatial_distribution
  if (is.null(df)) return(tibble())
  df %>%
    mutate(
      county_code = as.character(county_code),
      site_run = x$site_run,
      site_name = x$site_name,
      local_cluster = as.integer(traj_cluster_ra),
      county_fips = sprintf("%05d", readr::parse_number(county_code))
    )
})

county_cluster_pooled <- spatial_distribution %>%
  left_join(cluster_map %>% select(site_run, local_cluster, pooled_cluster, severity_score), by = c("site_run", "local_cluster")) %>%
  filter(!is.na(pooled_cluster)) %>%
  group_by(county_fips, pooled_cluster) %>%
  summarise(
    n = sum(n, na.rm = TRUE),
    severity_score = mean(severity_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(county_fips) %>%
  mutate(
    county_total = sum(n, na.rm = TRUE),
    prop = n / county_total
  ) %>%
  ungroup()

county_summary_pooled <- county_cluster_pooled %>%
  group_by(county_fips) %>%
  summarise(
    county_total = sum(n, na.rm = TRUE),
    dominant_cluster = pooled_cluster[which.max(n)],
    dominant_cluster_prop = max(prop, na.rm = TRUE),
    worst_cluster_prop = sum(n[pooled_cluster == max(pooled_cluster)], na.rm = TRUE) / sum(n, na.rm = TRUE),
    mean_severity_rank = weighted.mean(pooled_cluster, n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    prototype_key %>%
      transmute(
        pooled_cluster,
        pooled_cluster_label,
        dominant_cluster_phenotype = phenotype_label_short
      ),
    by = c("dominant_cluster" = "pooled_cluster")
  )

write_csv_safe(county_summary_pooled, file.path(OUT_DIR, "pooled_county_cluster_summary.csv"))

county_fips_lookup <- maps::county.fips %>%
  mutate(county_fips = sprintf("%05d", fips))

county_map_df <- ggplot2::map_data("county") %>%
  mutate(polyname = paste(region, subregion, sep = ",")) %>%
  left_join(county_fips_lookup, by = "polyname") %>%
  filter(!substr(county_fips, 1, 2) %in% c("02", "15", "72")) %>%
  left_join(county_summary_pooled, by = "county_fips")

cluster_palette <- c(
  "1" = "#d8e4bc",
  "2" = "#f4d06f",
  "3" = "#7fb7be",
  "4" = "#d96c75",
  "5" = "#7a0019"
)

phenotype_levels <- county_summary_pooled %>%
  distinct(dominant_cluster, dominant_cluster_phenotype) %>%
  arrange(dominant_cluster) %>%
  pull(dominant_cluster_phenotype)

cluster_palette_named <- stats::setNames(unname(cluster_palette[as.character(seq_along(phenotype_levels))]), phenotype_levels)

p_map_dom <- ggplot(
  county_map_df,
  aes(
    long,
    lat,
    group = group,
    fill = factor(dominant_cluster_phenotype, levels = phenotype_levels)
  )
) +
  geom_polygon(color = NA) +
  coord_quickmap() +
  scale_fill_manual(
    values = cluster_palette_named,
    na.value = "gray92",
    name = "Dominant pooled cluster phenotype",
    labels = \(x) stringr::str_wrap(x, width = 22)
  ) +
  labs(title = "Dominant Consensus Trajectory Cluster by County") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.key.width = unit(0.8, "cm"),
    legend.key.height = unit(1.1, "cm"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold")
  )

p_map_worst <- ggplot(county_map_df, aes(long, lat, group = group, fill = worst_cluster_prop)) +
  geom_polygon(color = NA) +
  coord_quickmap() +
  scale_fill_gradientn(
    colors = c("#f7fbff", "#9ecae1", "#3182bd", "#08519c"),
    na.value = "gray92",
    labels = scales::percent_format(accuracy = 1),
    name = "Share of county patients in worst pooled cluster",
    guide = guide_colorbar(
      title.position = "top",
      barheight = unit(10, "cm"),
      barwidth = unit(0.9, "cm")
    )
  ) +
  labs(title = "Share of County Patients in Worst Cluster") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

p_map_sev <- ggplot(county_map_df, aes(long, lat, group = group, fill = mean_severity_rank)) +
  geom_polygon(color = NA) +
  coord_quickmap() +
  scale_fill_gradientn(
    colors = c("#eef5db", "#b8d8ba", "#72b7b2", "#f28e2b", "#b22222"),
    na.value = "gray92",
    name = "Mean pooled trajectory severity rank",
    guide = guide_colorbar(
      title.position = "top",
      barheight = unit(10, "cm"),
      barwidth = unit(0.9, "cm")
    )
  ) +
  labs(title = "County-Level Mean Trajectory Severity") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

save_plot_safe(p_map_dom, file.path(OUT_DIR, "map_conus_dominant_pooled_cluster.png"), width = 13, height = 8)
save_plot_safe(p_map_worst, file.path(OUT_DIR, "map_conus_worst_cluster_share.png"), width = 13, height = 8)
save_plot_safe(p_map_sev, file.path(OUT_DIR, "map_conus_mean_cluster_severity.png"), width = 13, height = 8)

run_inventory <- tibble(
  site_run = map_chr(site_payload, "site_run"),
  site_name = map_chr(site_payload, "site_name"),
  has_centroids = map_lgl(site_payload, ~ !is.null(.x$centroids)),
  has_severity = map_lgl(site_payload, ~ !is.null(.x$severity)),
  has_static = map_lgl(site_payload, ~ !is.null(.x$static)),
  has_model_outputs = map_lgl(site_payload, ~ !is.null(.x$model_outputs)),
  has_overall_cont = map_lgl(site_payload, ~ !is.null(.x$overall_cont)),
  has_overall_cat = map_lgl(site_payload, ~ !is.null(.x$overall_cat)),
  has_bycluster_cont = map_lgl(site_payload, ~ !is.null(.x$bycluster_cont)),
  has_bycluster_cat = map_lgl(site_payload, ~ !is.null(.x$bycluster_cat)),
  has_contrast_no2 = map_lgl(site_payload, ~ !is.null(.x$contrast_no2)),
  has_contrast_pm25 = map_lgl(site_payload, ~ !is.null(.x$contrast_pm25))
)

write_csv_safe(run_inventory, file.path(OUT_DIR, "pooled_run_inventory.csv"))

notes <- tibble(
  item = c(
    "excluded_sites",
    "cluster_matching_method",
    "continuous_table1_pooling",
    "categorical_table1_pooling",
    "formal_meta_models",
    "cluster_specific_exposure_pooling",
    "known_input_gap",
    "interaction_forest_method",
    "map_method"
  ),
  detail = c(
    paste(EXCLUDED_SITES, collapse = "; "),
    "Site clusters were matched to 5 consensus prototypes using iterative minimum-distance assignment on exported centroid features with one-to-one matching within site.",
    "Continuous Table 1 summaries were pooled from site means, SDs, and sample sizes. Medians and quartiles are left blank because patient-level pooled quantiles are not available in the federated export.",
    "Categorical Table 1 summaries were pooled by summing site counts and denominators.",
    "Formal random-effects meta-analysis was limited to model terms that are label-stable across sites.",
    "Cluster-specific exposure summaries use the exported modeled low-vs-high contrast tables after consensus cluster remapping; these are weighted summaries rather than inverse-variance meta-analyses because standard errors were not exported for the contrasts.",
    "Michigan is missing cluster_signature_rs_hourly and cluster_signature_arf_hourly, but the files needed for the current pooling workflow were present.",
    "Cluster interaction forest plots and pooled interaction summaries use the site-level interaction coefficients from glm_interaction_no2 and glm_interaction_pm25. Terms were relabeled to the matching pooled trajectory phenotype using the consensus cluster mapping before random-effects pooling.",
    "CONUS maps aggregate county-level cluster distributions across sites after remapping local clusters to pooled consensus clusters."
  )
)

write_csv_safe(notes, file.path(OUT_DIR, "pooling_notes.csv"))

message("Pooled outputs written to: ", OUT_DIR)
