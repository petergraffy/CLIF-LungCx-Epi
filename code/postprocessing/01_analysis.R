suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

ROOT_DIR <- "."
SITE_DIR <- file.path(ROOT_DIR, "sites")
MODEL_EXCLUDED_SITES <- character(0)
RUN_DIRS_ALL <- sort(list.dirs(SITE_DIR, recursive = FALSE, full.names = TRUE))
RUN_DIRS_ALL <- RUN_DIRS_ALL[file.info(RUN_DIRS_ALL)$isdir %in% TRUE]
RUN_DATE <- format(Sys.Date(), "%Y%m%d")
OUT_DIR <- file.path(ROOT_DIR, "output", paste0("pooled_", RUN_DATE))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

write_csv_safe <- function(df, path) {
  readr::write_csv(df, path, na = "")
}

save_plot_safe <- function(plot_obj, path, width = 10, height = 7, dpi = 320) {
  ggplot2::ggsave(filename = path, plot = plot_obj, width = width, height = height, dpi = dpi, bg = "white")
}

save_diagram_safe <- function(diagram_obj, path, width = 12, height = 8.5) {
  svg_txt <- DiagrammeRsvg::export_svg(diagram_obj)
  rsvg::rsvg_png(
    charToRaw(svg_txt),
    file = path,
    width = width * 320,
    height = height * 320
  )
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
    "MINIMAL_ARF" = "Minimal ARF",
    "STABLE_ARF" = "Stable ARF",
    "RESOLVING_ARF" = "Resolving ARF",
    "PERSISTENT_ARF" = "Persistent ARF",
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
    "Minimal ARF" = "Minimal ARF",
    "Stable ARF" = "Stable ARF",
    "Resolving ARF" = "Resolving ARF",
    "Persistent ARF" = "Persistent ARF",
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

harmonize_model_term <- function(term) {
  term <- as.character(term)
  dplyr::case_when(
    term == "age_years" ~ "age_years",
    term == "admit_year" ~ "admit_year",
    term == "charlson_score" ~ "charlson_score",
    term == "advanced_cancer_any_poaTRUE" ~ "advanced_cancer_any_poaTRUE",
    stringr::str_to_lower(term) == "sex_categorymale" ~ "sex_categoryMale",
    stringr::str_to_lower(term) == "sex_categoryunknown" ~ "sex_categoryUnknown",
    stringr::str_detect(stringr::str_to_lower(term), "^race_category") ~ {
      lvl <- stringr::str_remove(term, "^race_category")
      lvl_low <- stringr::str_to_lower(lvl)
      paste0(
        "race_category",
        dplyr::case_when(
          lvl_low == "asian" ~ "Asian",
          lvl_low == "black or african american" ~ "Black or African American",
          lvl_low == "native hawaiian or other pacific islander" ~ "Native Hawaiian or Other Pacific Islander",
          lvl_low == "other" ~ "Other",
          lvl_low == "unknown" ~ "Unknown",
          lvl_low == "white" ~ "White",
          TRUE ~ lvl
        )
      )
    },
    TRUE ~ term
  )
}

pretty_covariate <- function(term) {
  dplyr::recode(
    term,
    "age_years" = "Age, per year",
    "admit_year" = "Admission year, per year",
    "charlson_score" = "Charlson score, per point",
    "advanced_cancer_any_poaTRUE" = "Advanced cancer present on admission",
    "sex_categoryMale" = "Male sex",
    "sex_categoryUnknown" = "Unknown sex",
    "race_categoryAsian" = "Asian race",
    "race_categoryBlack or African American" = "Black race",
    "race_categoryNative Hawaiian or Other Pacific Islander" = "Native Hawaiian or Pacific Islander race",
    "race_categoryOther" = "Other race",
    "race_categoryUnknown" = "Unknown race",
    "race_categoryWhite" = "White race",
    .default = term
  )
}

read_stub_csv <- function(run_dir, stub) {
  hits <- Sys.glob(file.path(run_dir, paste0(stub, "*.csv")))
  if (!length(hits)) return(NULL)
  out <- readr::read_csv(hits[[1]], show_col_types = FALSE, col_types = cols(.default = col_guess()))
  parse_count_value <- function(x) {
    x_chr <- trimws(as.character(x))
    parsed <- suppressWarnings(readr::parse_number(x_chr))
    dplyr::case_when(
      is.na(x_chr) | x_chr == "" ~ NA_real_,
      stringr::str_detect(x_chr, "^<") & !is.na(parsed) ~ parsed / 2,
      TRUE ~ parsed
    )
  }
  numeric_candidates <- c(
    "traj_cluster_ra", "severity_rank", "n_cluster", "n", "denom", "prop",
    "n_nonmissing", "mean", "sd", "median", "p25", "p75", "estimate", "std_error",
    "statistic", "p_value", "effect", "conf_low", "conf_high", "low", "high",
    "risk_diff", "risk_ratio", "events"
  )
  present_numeric <- intersect(numeric_candidates, names(out))
  if (length(present_numeric)) {
    count_like <- intersect(c("traj_cluster_ra", "severity_rank", "n_cluster", "n", "denom", "n_nonmissing", "events"), present_numeric)
    other_numeric <- setdiff(present_numeric, count_like)
    if (length(count_like)) {
      out <- out %>% mutate(across(all_of(count_like), parse_count_value))
    }
    if (length(other_numeric)) {
      out <- out %>% mutate(across(all_of(other_numeric), ~ suppressWarnings(readr::parse_number(as.character(.x)))))
    }
  }
  out
}

repair_categorical_summary <- function(df) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(df)
  required_cols <- c("variable", "level", "n")
  if (!all(required_cols %in% names(df))) return(df)
  df <- df %>%
    mutate(
      n = suppressWarnings(as.numeric(n)),
      denom = suppressWarnings(as.numeric(denom))
    ) %>%
    group_by(variable) %>%
    mutate(
      denom = ifelse(is.na(denom), sum(n, na.rm = TRUE), denom),
      prop = ifelse(!is.na(denom) & denom > 0, n / denom, as.numeric(NA))
    ) %>%
    ungroup()
  df
}

reconstruct_overall_cat_from_bycluster <- function(df, site_run, site_name) {
  if (is.null(df) || !is.data.frame(df) || !nrow(df)) return(tibble())
  required_cols <- c("variable", "level", "summary_type", "stratum", "n", "denom")
  if (!all(required_cols %in% names(df))) return(tibble())

  stratum_denoms <- df %>%
    mutate(
      n = suppressWarnings(as.numeric(n)),
      denom = suppressWarnings(as.numeric(denom))
    ) %>%
    filter(!is.na(denom)) %>%
    group_by(variable, stratum) %>%
    summarise(denom_stratum = max(denom, na.rm = TRUE), .groups = "drop")

  if (!nrow(stratum_denoms)) return(tibble())

  overall_denom <- stratum_denoms %>%
    group_by(variable) %>%
    summarise(denom = sum(denom_stratum, na.rm = TRUE), .groups = "drop")

  df %>%
    mutate(
      n = suppressWarnings(as.numeric(n)),
      level = harmonize_level(level, variable)
    ) %>%
    group_by(variable, level, summary_type) %>%
    summarise(n = sum(n, na.rm = TRUE), .groups = "drop") %>%
    left_join(overall_denom, by = "variable") %>%
    mutate(
      prop = ifelse(!is.na(denom) & denom > 0, n / denom, as.numeric(NA)),
      site_run = site_run,
      site_name = site_name,
      stratum = "Overall"
    ) %>%
    relocate(site_name, stratum, variable, level, summary_type)
}

rescale_count_columns <- function(df, ratio) {
  if (is.null(df) || !is.data.frame(df) || is.na(ratio) || ratio <= 0 || ratio == 1) return(df)
  count_cols <- intersect(c("n_cluster", "n", "denom", "n_nonmissing", "events"), names(df))
  if (!length(count_cols)) return(df)
  df <- df %>%
    mutate(across(all_of(count_cols), ~ round(.x * ratio)))
  if (all(c("n", "denom") %in% names(df))) {
    df <- df %>%
      mutate(prop = dplyr::if_else(!is.na(denom) & denom > 0, n / denom, as.numeric(NA)))
  }
  df
}

infer_site_run <- function(run_dir) {
  basename(run_dir)
}

infer_site_name <- function(run_dir) {
  nm <- sub("^run_", "", basename(run_dir))
  sub("_[0-9]{8}$", "", nm)
}

RUN_DIRS <- RUN_DIRS_ALL
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
  cluster_gcs_summary = "cluster_gcs_summary",
  consort_flow = "consort_flow_counts",
  consort_reasons = "consort_exclusion_reasons",
  mediation_no2 = "mediation_style_decomposition_no2",
  mediation_no2_cluster = "mediation_style_decomposition_no2_cluster",
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
  if (!is.null(tables$overall_cat) && nrow(tables$overall_cat)) {
    tables$overall_cat <- repair_categorical_summary(tables$overall_cat)
  }
  if ((is.null(tables$overall_cat) || !nrow(tables$overall_cat)) && !is.null(tables$bycluster_cat) && nrow(tables$bycluster_cat)) {
    tables$overall_cat <- reconstruct_overall_cat_from_bycluster(tables$bycluster_cat, site_run, site_name)
  }
  tables$site_run <- site_run
  tables$site_name <- site_name
  tables
})
names(site_payload) <- map_chr(site_payload, "site_run")

hopkins_run <- names(site_payload)[map_chr(site_payload, "site_name") == "Hopkins"]
if (length(hopkins_run) == 1) {
  hopkins_flow <- site_payload[[hopkins_run]]$consort_flow
  hopkins_overall_cont <- site_payload[[hopkins_run]]$overall_cont
  hopkins_target_n <- hopkins_flow %>%
    filter(step == "Any ICU segment present") %>%
    summarise(n = first(n_remaining)) %>%
    pull(n)
  hopkins_source_n <- hopkins_overall_cont %>%
    filter(variable == "age_years", stratum == "Overall") %>%
    summarise(n = first(n_nonmissing)) %>%
    pull(n)
  hopkins_ratio <- hopkins_target_n / hopkins_source_n

  if (is.finite(hopkins_ratio) && hopkins_ratio > 0 && hopkins_ratio < 1) {
    rescaled_tables <- c(
      "centroids", "static", "cluster_gcs_summary", "rs_signature", "arf_signature",
      "cluster_exposure_mortality", "spatial_summary", "spatial_distribution",
      "overall_cont", "overall_cat", "bycluster_cont", "bycluster_cat"
    )
    for (nm in rescaled_tables) {
      site_payload[[hopkins_run]][[nm]] <- rescale_count_columns(site_payload[[hopkins_run]][[nm]], hopkins_ratio)
    }
  }
}

model_site_runs <- map_chr(site_payload, "site_run")[!map_chr(site_payload, "site_name") %in% MODEL_EXCLUDED_SITES]
site_payload_model <- site_payload[model_site_runs]

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
pooled_k_target <- centroids %>%
  count(site_run, name = "k_local") %>%
  summarise(k = round(median(k_local, na.rm = TRUE))) %>%
  pull(k)
pooled_k_target <- max(2L, min(as.integer(pooled_k_target %||% 5L), nrow(feature_frame)))
set.seed(20260408)
init_fit <- kmeans(feature_matrix, centers = pooled_k_target, nstart = 100)
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
    dominant_arf = c("NO_ARF", "HYPOX", "HYPER", "MIXED")[max.col(select(., pct_no_arf, pct_hypox, pct_hyper, pct_mixed), ties.method = "first")],
    arf_burden = 1 - pct_no_arf
  ) %>%
  arrange(severity_score) %>%
  mutate(
    pooled_cluster = row_number(),
    pooled_severity_rank = row_number(),
    is_worst_cluster = pooled_cluster == max(pooled_cluster),
    arf_label = case_when(
      dominant_rs == "IMV" | arf_burden >= 0.15 ~ "PERSISTENT_ARF",
      dominant_rs == "NIV" | arf_burden >= 0.04 ~ "STABLE_ARF",
      dominant_rs == "LOW_O2" ~ "RESOLVING_ARF",
      TRUE ~ "MINIMAL_ARF"
    ),
    pooled_cluster_label = paste0("Cluster ", pooled_cluster, ": ", label_arf(arf_label)),
    phenotype_label = label_arf(arf_label),
    phenotype_label_short = label_arf_short(label_arf(arf_label))
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

summary_anova_p <- function(df) {
  df <- df %>%
    filter(!is.na(n_nonmissing), !is.na(mean), !is.na(sd), n_nonmissing > 1)
  k <- nrow(df)
  if (k < 2) return(NA_real_)
  n_total <- sum(df$n_nonmissing)
  if (n_total <= k) return(NA_real_)
  grand_mean <- weighted.mean(df$mean, df$n_nonmissing)
  ss_within <- sum((df$n_nonmissing - 1) * (df$sd^2))
  ss_between <- sum(df$n_nonmissing * (df$mean - grand_mean)^2)
  df_between <- k - 1
  df_within <- n_total - k
  if (df_between <= 0 || df_within <= 0 || ss_within < 0) return(NA_real_)
  ms_between <- ss_between / df_between
  ms_within <- ss_within / df_within
  if (ms_within <= 0) return(NA_real_)
  stats::pf(ms_between / ms_within, df1 = df_between, df2 = df_within, lower.tail = FALSE)
}

categorical_prop_p <- function(df) {
  df <- df %>%
    filter(!is.na(n), !is.na(denom), denom >= n)
  if (nrow(df) < 2) return(NA_real_)
  mat <- cbind(selected = df$n, other = df$denom - df$n)
  if (any(rowSums(mat) == 0)) return(NA_real_)
  suppressWarnings(stats::chisq.test(mat)$p.value)
}

fmt_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

cluster_ns <- cluster_map %>%
  group_by(pooled_cluster, pooled_cluster_label) %>%
  summarise(n = sum(n, na.rm = TRUE), .groups = "drop")

overall_n <- pooled_overall_cont %>%
  filter(stratum == "Overall", variable == "age_years") %>%
  summarise(n = sum(n_nonmissing, na.rm = TRUE)) %>%
  pull(n)

cont_labels <- c(
  age_years = "Age, mean (SD)",
  charlson_score = "Charlson score, mean (SD)",
  sofa_total = "SOFA total, mean (SD)",
  pm25_5y = "PM₂.₅ 5-year exposure, mean (SD)",
  no2_5y = "NO₂ 5-year exposure, mean (SD)",
  icu_los_hours = "ICU LOS (hours), mean (SD)",
  imv_hours_72h = "IMV hours in first 72h, mean (SD)",
  vaso_hours_72h = "Vasopressor hours in first 72h, mean (SD)"
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
  "any_imv_72h", "1", "Any IMV in first 72h, n (%)",
  "vaso_any_72h", "1", "Any vasopressor in first 72h, n (%)"
)

site_order <- c(sort(unique(overall_cont$site_name)), "Overall")

site_table1_cont <- bind_rows(
  overall_cont %>%
    filter(variable %in% names(cont_labels), stratum == "Overall") %>%
    transmute(site_col = site_name, label = unname(cont_labels[variable]), value = fmt_mean_sd(mean, sd)),
  pooled_overall_cont %>%
    filter(variable %in% names(cont_labels), stratum == "Overall") %>%
    transmute(site_col = "Overall", label = unname(cont_labels[variable]), value = fmt_mean_sd(mean, sd))
) %>%
  mutate(site_col = factor(site_col, levels = site_order)) %>%
  arrange(site_col) %>%
  pivot_wider(names_from = site_col, values_from = value)

site_table1_cat <- bind_rows(
  overall_cat %>%
    inner_join(cat_specs, by = c("variable", "level")) %>%
    transmute(site_col = site_name, label, value = fmt_n_pct(n, denom)),
  pooled_overall_cat %>%
    inner_join(cat_specs, by = c("variable", "level")) %>%
    transmute(site_col = "Overall", label, value = fmt_n_pct(n, denom))
) %>%
  mutate(site_col = factor(site_col, levels = site_order)) %>%
  arrange(site_col) %>%
  pivot_wider(names_from = site_col, values_from = value)

site_table1_n <- bind_rows(
  overall_cont %>%
    filter(variable == "age_years", stratum == "Overall") %>%
    transmute(site_col = site_name, label = "N", value = scales::comma(n_nonmissing)),
  tibble(site_col = "Overall", label = "N", value = scales::comma(overall_n))
) %>%
  mutate(site_col = factor(site_col, levels = site_order)) %>%
  arrange(site_col) %>%
  pivot_wider(names_from = site_col, values_from = value)

manuscript_table1 <- bind_rows(
  site_table1_n,
  site_table1_cont,
  site_table1_cat
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

table1_cat_formatted_specs <- tribble(
  ~Characteristic, ~Level, ~variable, ~level_value,
  "Sex", "Female", "sex_category", "Female",
  "Sex", "Male", "sex_category", "Male",
  "Race", "American Indian or Alaska Native", "race_category", "American Indian or Alaska Native",
  "Race", "Asian", "race_category", "Asian",
  "Race", "Black or African American", "race_category", "Black or African American",
  "Race", "Native Hawaiian or Other Pacific Islander", "race_category", "Native Hawaiian or Other Pacific Islander",
  "Race", "Other", "race_category", "Other",
  "Race", "White", "race_category", "White",
  "Ethnicity", "Hispanic", "ethnicity_category", "Hispanic",
  "Ethnicity", "Non-Hispanic", "ethnicity_category", "Non-Hispanic",
  "Advanced cancer present on admission", "No", "advanced_cancer_any_poa", "0",
  "Advanced cancer present on admission", "Yes", "advanced_cancer_any_poa", "1",
  "In-hospital death", "No", "death_in_hosp", "0",
  "In-hospital death", "Yes", "death_in_hosp", "1",
  "Hospice discharge", "No", "hospice_discharge", "0",
  "Hospice discharge", "Yes", "hospice_discharge", "1",
  "Death or hospice discharge", "No", "death_or_hospice", "0",
  "Death or hospice discharge", "Yes", "death_or_hospice", "1",
  "Any invasive mechanical ventilation in first 72 hours", "No", "any_imv_72h", "0",
  "Any invasive mechanical ventilation in first 72 hours", "Yes", "any_imv_72h", "1",
  "Any vasopressor in first 72 hours", "No", "vaso_any_72h", "0",
  "Any vasopressor in first 72 hours", "Yes", "vaso_any_72h", "1"
)

table1_formatted_n <- bind_rows(
  overall_cont %>%
    filter(variable == "age_years", stratum == "Overall") %>%
    transmute(site_col = site_name, Characteristic = "N", Level = "", value = scales::comma(n_nonmissing)),
  tibble(site_col = "Overall", Characteristic = "N", Level = "", value = scales::comma(overall_n))
) %>%
  mutate(site_col = factor(site_col, levels = site_order)) %>%
  arrange(site_col) %>%
  pivot_wider(names_from = site_col, values_from = value)

table1_formatted_cont <- bind_rows(
  overall_cont %>%
    filter(variable %in% names(cont_labels), stratum == "Overall") %>%
    transmute(site_col = site_name, Characteristic = unname(cont_labels[variable]), Level = "", value = fmt_mean_sd(mean, sd)),
  pooled_overall_cont %>%
    filter(variable %in% names(cont_labels), stratum == "Overall") %>%
    transmute(site_col = "Overall", Characteristic = unname(cont_labels[variable]), Level = "", value = fmt_mean_sd(mean, sd))
) %>%
  mutate(site_col = factor(site_col, levels = site_order)) %>%
  arrange(site_col) %>%
  pivot_wider(names_from = site_col, values_from = value)

table1_formatted_cat <- bind_rows(
  overall_cat %>%
    inner_join(table1_cat_formatted_specs, by = c("variable", "level" = "level_value")) %>%
    transmute(site_col = site_name, Characteristic, Level, value = fmt_n_pct(n, denom)),
  pooled_overall_cat %>%
    inner_join(table1_cat_formatted_specs, by = c("variable", "level" = "level_value")) %>%
    transmute(site_col = "Overall", Characteristic, Level, value = fmt_n_pct(n, denom))
) %>%
  mutate(site_col = factor(site_col, levels = site_order)) %>%
  arrange(site_col) %>%
  pivot_wider(names_from = site_col, values_from = value)

table1_formatted_order <- c(
  "N",
  unname(cont_labels),
  "Sex",
  "Race",
  "Ethnicity",
  "Advanced cancer present on admission",
  "In-hospital death",
  "Hospice discharge",
  "Death or hospice discharge",
  "Any invasive mechanical ventilation in first 72 hours",
  "Any vasopressor in first 72 hours"
)

manuscript_table1_formatted <- bind_rows(
  table1_formatted_n,
  table1_formatted_cont,
  table1_formatted_cat
) %>%
  mutate(
    sort_order = match(Characteristic, table1_formatted_order),
    level_order = case_when(
      Level == "" ~ 0,
      Level %in% c("Female", "Hispanic", "No") ~ 1,
      Level %in% c("Male", "Non-Hispanic", "Yes") ~ 2,
      TRUE ~ 3
    )
  ) %>%
  arrange(sort_order, level_order, Level) %>%
  select(-sort_order, -level_order)

write_csv_safe(manuscript_table1_formatted, file.path(OUT_DIR, "manuscript_table1_formatted.csv"))

manuscript_table1_overall_formatted <- manuscript_table1_formatted %>%
  select(Characteristic, Level, Overall)

write_csv_safe(manuscript_table1_overall_formatted, file.path(OUT_DIR, "manuscript_table1_overall_formatted.csv"))

table1_cols <- names(manuscript_table1)
blank_row <- as.list(rep(NA_character_, length(table1_cols)))
names(blank_row) <- table1_cols

table1_section_cont <- blank_row
table1_section_cont$label <- "Continuous variables"

table1_section_cat <- blank_row
table1_section_cat$label <- "Categorical variables"

manuscript_table1_polished <- bind_rows(
  manuscript_table1 %>% filter(label == "N"),
  tibble::as_tibble_row(table1_section_cont),
  manuscript_table1 %>% filter(label %in% unname(cont_labels)),
  tibble::as_tibble_row(table1_section_cat),
  manuscript_table1 %>% filter(label %in% cat_specs$label)
)

write_csv_safe(manuscript_table1_polished, file.path(OUT_DIR, "manuscript_table1_polished.csv"))

table2_cont_p <- pooled_bycluster_cont %>%
  filter(variable %in% names(cont_labels)) %>%
  group_by(variable) %>%
  summarise(p_value = summary_anova_p(pick(everything())), .groups = "drop")

table2_cat_p <- pooled_bycluster_cat %>%
  inner_join(cat_specs, by = c("variable", "level")) %>%
  group_by(label) %>%
  summarise(p_value = categorical_prop_p(pick(everything())), .groups = "drop")

cluster_col_order <- paste0("Cluster ", sort(unique(cluster_ns$pooled_cluster)))

table2_n <- cluster_ns %>%
  transmute(
    label = "N",
    pooled_cluster,
    value = scales::comma(n)
  ) %>%
  pivot_wider(names_from = pooled_cluster, values_from = value, names_prefix = "Cluster ") %>%
  {
    missing_cols <- setdiff(cluster_col_order, names(.))
    if (length(missing_cols)) mutate(., across(any_of(names(.)), identity), !!!setNames(rep(list(NA_character_), length(missing_cols)), missing_cols)) else .
  } %>%
  select(label, all_of(cluster_col_order)) %>%
  mutate(`P value` = NA_character_)

table2_cont <- pooled_bycluster_cont %>%
  filter(variable %in% names(cont_labels)) %>%
  transmute(
    variable,
    pooled_cluster,
    label = unname(cont_labels[variable]),
    value = fmt_mean_sd(mean, sd)
  ) %>%
  pivot_wider(names_from = pooled_cluster, values_from = value, names_prefix = "Cluster ") %>%
  left_join(table2_cont_p, by = "variable") %>%
  {
    missing_cols <- setdiff(cluster_col_order, names(.))
    if (length(missing_cols)) mutate(., !!!setNames(rep(list(NA_character_), length(missing_cols)), missing_cols)) else .
  } %>%
  transmute(label, !!!rlang::syms(cluster_col_order), `P value` = fmt_p(p_value))

table2_cat <- pooled_bycluster_cat %>%
  inner_join(cat_specs, by = c("variable", "level")) %>%
  transmute(
    label,
    pooled_cluster,
    value = fmt_n_pct(n, denom)
  ) %>%
  pivot_wider(names_from = pooled_cluster, values_from = value, names_prefix = "Cluster ") %>%
  left_join(table2_cat_p, by = "label") %>%
  {
    missing_cols <- setdiff(cluster_col_order, names(.))
    if (length(missing_cols)) mutate(., !!!setNames(rep(list(NA_character_), length(missing_cols)), missing_cols)) else .
  } %>%
  transmute(label, !!!rlang::syms(cluster_col_order), `P value` = fmt_p(p_value))

table2_order <- c(
  "N",
  unname(cont_labels),
  cat_specs$label
)

manuscript_table2 <- bind_rows(
  table2_n,
  table2_cont,
  table2_cat
) %>%
  mutate(sort_order = match(label, table2_order)) %>%
  arrange(sort_order) %>%
  select(-sort_order)

write_csv_safe(manuscript_table2, file.path(OUT_DIR, "manuscript_table2_by_cluster.csv"))

model_outputs <- map_dfr(site_payload_model, function(x) {
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

covariate_models <- tribble(
  ~model, ~term_group,
  "glm_mortality", "Mortality",
  "glm_los_quasipoisson", "ICU LOS",
  "cox72_landmark", "Post-72h survival",
  "polr_severity_rank", "Trajectory severity"
)

candidate_covariate_terms <- c(
  "age_years",
  "admit_year",
  "charlson_score",
  "advanced_cancer_any_poaTRUE",
  "sex_categoryMale",
  "sex_categoryUnknown",
  "race_categoryAsian",
  "race_categoryBlack or African American",
  "race_categoryNative Hawaiian or Other Pacific Islander",
  "race_categoryOther",
  "race_categoryUnknown",
  "race_categoryWhite"
)

outcome_covariates <- c(
  "Age, per year",
  "Charlson score, per point",
  "Advanced cancer present on admission",
  "Male sex"
)

race_covariates <- c(
  "Asian race",
  "Black race",
  "White race",
  "Other race"
)

pooled_covariate_meta <- model_outputs %>%
  mutate(term_clean = harmonize_model_term(term)) %>%
  inner_join(covariate_models, by = "model") %>%
  filter(term_clean %in% candidate_covariate_terms) %>%
  group_by(term_group, model, term_clean, effect_type) %>%
  group_modify(~ {
    effect_type_val <- .y$effect_type[[1]] %||% NA_character_
    meta_random_effects(.x) %>%
      mutate(
        pooled_effect_random = ifelse(!is.na(effect_type_val) && effect_type_val != "", exp(estimate_random), NA_real_),
        pooled_conf_low_random = ifelse(!is.na(effect_type_val) && effect_type_val != "", exp(conf_low_random), NA_real_),
        pooled_conf_high_random = ifelse(!is.na(effect_type_val) && effect_type_val != "", exp(conf_high_random), NA_real_),
        covariate_label = pretty_covariate(.y$term_clean[[1]]),
        site_included = paste(sort(unique(.x$site_name)), collapse = "; ")
      )
  }) %>%
  ungroup() %>%
  arrange(term_group, covariate_label)

write_csv_safe(pooled_covariate_meta, file.path(OUT_DIR, "pooled_covariate_meta_analysis.csv"))

race_covariate_table <- pooled_covariate_meta %>%
  filter(covariate_label %in% race_covariates) %>%
  transmute(
    outcome = term_group,
    race_covariate = covariate_label,
    pooled_effect = pooled_effect_random,
    conf_low = pooled_conf_low_random,
    conf_high = pooled_conf_high_random,
    estimate_ci = fmt_effect(pooled_effect_random, pooled_conf_low_random, pooled_conf_high_random),
    k,
    site_included
  ) %>%
  arrange(race_covariate, outcome)

write_csv_safe(race_covariate_table, file.path(OUT_DIR, "pooled_race_covariate_effects.csv"))

cluster_membership_covariates <- model_outputs %>%
  filter(model == "multinom_traj_cluster") %>%
  mutate(
    term_clean = harmonize_model_term(term),
    local_cluster = suppressWarnings(as.integer(outcome_level))
  ) %>%
  filter(term_clean %in% candidate_covariate_terms) %>%
  left_join(
    cluster_map %>% select(site_run, local_cluster, pooled_cluster, pooled_cluster_label),
    by = c("site_run", "local_cluster")
  ) %>%
  filter(!is.na(pooled_cluster)) %>%
  group_by(pooled_cluster, pooled_cluster_label, term_clean, effect_type) %>%
  group_modify(~ {
    effect_type_val <- .y$effect_type[[1]] %||% NA_character_
    meta_random_effects(.x) %>%
      mutate(
        pooled_effect_random = ifelse(!is.na(effect_type_val) && effect_type_val != "", exp(estimate_random), NA_real_),
        pooled_conf_low_random = ifelse(!is.na(effect_type_val) && effect_type_val != "", exp(conf_low_random), NA_real_),
        pooled_conf_high_random = ifelse(!is.na(effect_type_val) && effect_type_val != "", exp(conf_high_random), NA_real_),
        covariate_label = pretty_covariate(.y$term_clean[[1]]),
        site_included = paste(sort(unique(.x$site_name)), collapse = "; ")
      )
  }) %>%
  ungroup() %>%
  arrange(pooled_cluster, covariate_label)

write_csv_safe(cluster_membership_covariates, file.path(OUT_DIR, "pooled_cluster_membership_covariates.csv"))

mediation_no2 <- map_dfr(site_payload_model, function(x) {
  df <- x$mediation_no2
  if (is.null(df)) return(tibble())
  df %>% mutate(site_run = x$site_run, site_name = x$site_name, mediator = "Severity rank")
})

mediation_no2_cluster <- map_dfr(site_payload_model, function(x) {
  df <- x$mediation_no2_cluster
  if (is.null(df)) return(tibble())
  df %>% mutate(site_run = x$site_run, site_name = x$site_name, mediator = "Trajectory phenotype")
})

mediation_site_level <- bind_rows(mediation_no2, mediation_no2_cluster) %>%
  mutate(site_name = as.character(site_name))

mediation_pooled_models <- pooled_meta %>%
  filter(model %in% c("glm_total_no2", "glm_direct_no2_severity", "glm_direct_no2_cluster"), term == "no2_5y_z") %>%
  transmute(
    model,
    mediator = case_when(
      model == "glm_total_no2" ~ "Total effect",
      model == "glm_direct_no2_severity" ~ "Severity rank",
      model == "glm_direct_no2_cluster" ~ "Trajectory phenotype"
    ),
    pooled_beta = estimate_random,
    pooled_beta_low = conf_low_random,
    pooled_beta_high = conf_high_random,
    pooled_or = pooled_effect_random,
    pooled_or_low = pooled_conf_low_random,
    pooled_or_high = pooled_conf_high_random
  )

mediation_summary <- mediation_site_level %>%
  group_by(mediator) %>%
  summarise(
    n_sites = n_distinct(site_name),
    beta_total_mean = mean(beta_total, na.rm = TRUE),
    beta_direct_mean = mean(beta_direct, na.rm = TRUE),
    beta_indirect_mean = mean(beta_indirect_approx, na.rm = TRUE),
    pct_mediated_median = median(pct_mediated_approx, na.rm = TRUE),
    pct_mediated_iqr_low = quantile(pct_mediated_approx, 0.25, na.rm = TRUE),
    pct_mediated_iqr_high = quantile(pct_mediated_approx, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(mediation_pooled_models, by = "mediator")

write_csv_safe(mediation_site_level, file.path(OUT_DIR, "publication_table_mediation_site_level.csv"))
write_csv_safe(mediation_summary, file.path(OUT_DIR, "publication_table_mediation_summary.csv"))

mediation_plot_df <- mediation_site_level %>%
  transmute(
    site_name,
    mediator,
    `Total effect` = beta_total,
    `Direct effect` = beta_direct,
    `Indirect effect` = beta_indirect_approx
  ) %>%
  pivot_longer(cols = c(`Total effect`, `Direct effect`, `Indirect effect`), names_to = "effect_type", values_to = "beta") %>%
  mutate(effect_type = factor(effect_type, levels = c("Total effect", "Direct effect", "Indirect effect")))

p_mediation_site <- ggplot(mediation_plot_df, aes(x = effect_type, y = beta, group = site_name, color = site_name)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray70") +
  geom_line(alpha = 0.30, linewidth = 0.7) +
  geom_point(size = 2.2) +
  facet_wrap(~ mediator) +
  labs(
    title = "Site-Level Approximate NO₂ Mediation Decomposition",
    x = NULL,
    y = "Log-odds coefficient for NO₂",
    color = "Site"
  ) +
  theme_manuscript() +
  theme(legend.position = "bottom")

save_plot_safe(p_mediation_site, file.path(OUT_DIR, "figure_mediation_site_level_no2.png"), width = 10.5, height = 5.8)

total_beta_row <- mediation_pooled_models %>% filter(mediator == "Total effect") %>% slice(1)

mediation_pooled_plot_df <- bind_rows(
  mediation_pooled_models %>%
    filter(mediator != "Total effect") %>%
    transmute(
      mediator,
      effect_type = "Total effect",
      estimate = total_beta_row$pooled_beta[[1]],
      conf_low = total_beta_row$pooled_beta_low[[1]],
      conf_high = total_beta_row$pooled_beta_high[[1]]
    ),
  mediation_pooled_models %>%
    filter(mediator != "Total effect") %>%
    transmute(
      mediator,
      effect_type = "Direct effect",
      estimate = pooled_beta,
      conf_low = pooled_beta_low,
      conf_high = pooled_beta_high
    ),
  mediation_summary %>%
    transmute(
      mediator,
      effect_type = "Indirect effect",
      estimate = beta_indirect_mean,
      conf_low = NA_real_,
      conf_high = NA_real_
    )
) %>%
  mutate(effect_type = factor(effect_type, levels = c("Total effect", "Direct effect", "Indirect effect")))

pct_label_df <- mediation_summary %>%
  transmute(
    mediator,
    pct_label = sprintf(
      "Median mediated: %.0f%% (IQR %.0f%% to %.0f%%)",
      100 * pct_mediated_median,
      100 * pct_mediated_iqr_low,
      100 * pct_mediated_iqr_high
    )
  )

p_mediation_pooled <- ggplot(mediation_pooled_plot_df, aes(x = effect_type, y = estimate, color = mediator, group = mediator)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray70") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3.2) +
  geom_errorbar(
    data = ~ dplyr::filter(.x, !is.na(conf_low), !is.na(conf_high)),
    aes(ymin = conf_low, ymax = conf_high),
    width = 0.12,
    linewidth = 0.8
  ) +
  labs(
    title = "Pooled Attenuation of the NO₂ Association After Mediator Adjustment",
    x = NULL,
    y = "NO₂ log-odds coefficient",
    color = "Mediator"
  ) +
  theme_manuscript() +
  theme(legend.position = "bottom")

save_plot_safe(p_mediation_pooled, file.path(OUT_DIR, "figure_mediation_pooled_no2.png"), width = 9.5, height = 5.5)

p_mediation_combo <- p_mediation_site / p_mediation_pooled +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1.5, 1))

save_plot_safe(p_mediation_combo, file.path(OUT_DIR, "figure_mediation_no2_two_panel.png"), width = 11, height = 11)

contrast_no2 <- map_dfr(site_payload_model, function(x) {
  df <- x$contrast_no2
  if (is.null(df)) return(tibble())
  df %>%
    mutate(
      site_run = x$site_run,
      site_name = x$site_name,
      local_cluster = as.integer(traj_cluster_ra)
    )
})

contrast_pm25 <- map_dfr(site_payload_model, function(x) {
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
    outcome = case_when(
      term_group == "mortality" ~ "Mortality",
      term_group == "icu_los" ~ "ICU LOS",
      term_group == "post72_survival" ~ "Post-72h survival",
      term_group == "trajectory_severity" ~ "Trajectory severity",
      TRUE ~ NA_character_
    ),
    effect_display = fmt_effect(pooled_effect_random, pooled_conf_low_random, pooled_conf_high_random)
  ) %>%
  filter(!is.na(pooled_effect_random), !is.na(outcome)) %>%
  mutate(
    outcome = factor(outcome, levels = rev(c("ICU LOS", "Mortality", "Post-72h survival", "Trajectory severity"))),
    exposure = factor(exposure, levels = c("NO₂", "PM₂.₅")),
    y_base = c(
      "ICU LOS" = 4,
      "Mortality" = 3,
      "Post-72h survival" = 2,
      "Trajectory severity" = 1
    )[as.character(outcome)],
    y_pos = y_base + ifelse(exposure == "NO₂", 0.12, -0.12),
    label_x = pooled_conf_high_random * 1.03
  )

meta_xmin <- min(meta_plot_df$pooled_conf_low_random, na.rm = TRUE) * 0.97
meta_xmax <- max(meta_plot_df$label_x, na.rm = TRUE) * 1.04
meta_xmin <- min(meta_xmin, 0.8)
meta_xmax <- max(meta_xmax, 1.4)

p_meta <- ggplot(
  meta_plot_df,
  aes(x = pooled_effect_random, y = y_pos, xmin = pooled_conf_low_random, xmax = pooled_conf_high_random, color = exposure)
) +
  geom_vline(xintercept = 1, linetype = 2, color = "gray45") +
  geom_errorbar(width = 0.10, linewidth = 0.8, orientation = "y") +
  geom_point(size = 3.4) +
  geom_text(aes(x = label_x, label = effect_display), hjust = 0, size = 3.8, color = "gray15", show.legend = FALSE) +
  scale_x_log10(
    breaks = c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4),
    labels = c("0.8", "0.9", "1.0", "1.1", "1.2", "1.3", "1.4")
  ) +
  scale_y_continuous(
    breaks = c(4, 3, 2, 1),
    labels = c("ICU LOS", "Mortality", "Post-72h survival", "Trajectory severity"),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  labs(
    title = "Pooled Exposure Associations Across Outcomes",
    x = "Pooled effect estimate (random-effects, log scale)",
    y = NULL,
    color = NULL
  ) +
  coord_cartesian(xlim = c(meta_xmin, meta_xmax), clip = "off") +
  theme_manuscript() +
  theme(
    plot.margin = margin(12, 130, 12, 12),
    legend.position = "bottom",
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12)
  )

save_plot_safe(p_meta, file.path(OUT_DIR, "forest_overall_meta_effects.png"), width = 12.5, height = 7.5)

cluster_covariates <- c(
  "Age, per year",
  "Charlson score, per point",
  "Advanced cancer present on admission",
  "Male sex"
)

covariate_outcome_plot_df <- pooled_covariate_meta %>%
  filter(covariate_label %in% outcome_covariates) %>%
  mutate(
    covariate_label = factor(covariate_label, levels = rev(outcome_covariates)),
    term_group = factor(term_group, levels = c("ICU LOS", "Mortality", "Post-72h survival", "Trajectory severity")),
    y_base = c(
      "Age, per year" = 4,
      "Charlson score, per point" = 3,
      "Advanced cancer present on admission" = 2,
      "Male sex" = 1
    )[as.character(covariate_label)],
    y_pos = y_base + c(
      "ICU LOS" = 0.24,
      "Mortality" = 0.08,
      "Post-72h survival" = -0.08,
      "Trajectory severity" = -0.24
    )[as.character(term_group)],
    label_text = fmt_effect(pooled_effect_random, pooled_conf_low_random, pooled_conf_high_random),
    label_x = pooled_conf_high_random * 1.03
  )

covariate_outcome_xmax <- max(covariate_outcome_plot_df$pooled_conf_high_random, na.rm = TRUE) * 1.08
covariate_outcome_xmax <- max(covariate_outcome_xmax, 2.55)
covariate_outcome_xmax <- ceiling(covariate_outcome_xmax * 10) / 10

p_covariate_outcomes <- ggplot(
  covariate_outcome_plot_df,
  aes(x = pooled_effect_random, y = y_pos, xmin = pooled_conf_low_random, xmax = pooled_conf_high_random, color = term_group)
) +
  geom_vline(xintercept = 1, linetype = 2, color = "gray45") +
  geom_errorbar(width = 0.10, linewidth = 0.8, orientation = "y") +
  geom_point(size = 2.5) +
  geom_text(aes(x = label_x, label = label_text), hjust = 0, size = 2.8, color = "gray15", show.legend = FALSE) +
  scale_x_continuous(
    breaks = c(0.6, 0.8, 1.0, 1.1, 1.2, 1.5, 2.0, 2.5, 3.0),
    labels = c("0.6", "0.8", "1.0", "1.1", "1.2", "1.5", "2.0", "2.5", "3.0"),
    minor_breaks = c(0.9, 1.05, 1.15),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_y_continuous(
    breaks = c(4, 3, 2, 1),
    labels = outcome_covariates,
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_color_manual(
    values = c(
      "ICU LOS" = "#355070",
      "Mortality" = "#b56576",
      "Post-72h survival" = "#6d597a",
      "Trajectory severity" = "#2a9d8f"
    ),
    name = NULL
  ) +
  labs(
    title = "Pooled Covariate Associations Across Main Outcomes",
    x = "Pooled effect estimate (random-effects)",
    y = NULL
  ) +
  coord_cartesian(xlim = c(0.6, covariate_outcome_xmax), clip = "off") +
  theme_manuscript() +
  theme(
    plot.margin = margin(10, 95, 10, 10),
    legend.position = "bottom",
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(size = 10)
  )

save_plot_safe(p_covariate_outcomes, file.path(OUT_DIR, "forest_covariate_effects_outcomes.png"), width = 12, height = 8.5)

cluster_membership_plot_df <- cluster_membership_covariates %>%
  filter(covariate_label %in% cluster_covariates) %>%
  mutate(
    covariate_label = factor(covariate_label, levels = rev(cluster_covariates)),
    pooled_cluster_label = factor(
      pooled_cluster_label,
      levels = unique(prototype_key$pooled_cluster_label[order(prototype_key$pooled_cluster)])
    ),
    y_base = c(
      "Age, per year" = 4,
      "Charlson score, per point" = 3,
      "Advanced cancer present on admission" = 2,
      "Male sex" = 1
    )[as.character(covariate_label)],
    y_pos = y_base + c(-0.28, -0.14, 0, 0.14, 0.28)[as.integer(pooled_cluster_label)],
    label_text = fmt_effect(pooled_effect_random, pooled_conf_low_random, pooled_conf_high_random),
    label_x = pooled_conf_high_random * 1.03
  )

phenotype_palette_master <- c(
  "Minimal ARF" = "#2E8B57",
  "Stable ARF" = "#D4A017",
  "Resolving ARF" = "#2B6CB0",
  "Persistent ARF" = "#C05621"
)

cluster_plot_palette <- prototype_key %>%
  distinct(pooled_cluster, pooled_cluster_label, phenotype_label_short) %>%
  arrange(pooled_cluster) %>%
  mutate(color = phenotype_palette_master[phenotype_label_short]) %>%
  { stats::setNames(.$color, .$pooled_cluster_label) }

cluster_plot_palette_short <- prototype_key %>%
  distinct(pooled_cluster, phenotype_label_short) %>%
  arrange(pooled_cluster) %>%
  mutate(
    cluster_label = paste0("C", pooled_cluster, ": ", phenotype_label_short),
    color = phenotype_palette_master[phenotype_label_short]
  ) %>%
  { stats::setNames(.$color, .$cluster_label) }

cluster_plot_palette_phenotype_short <- prototype_key %>%
  distinct(pooled_cluster, phenotype_label_short) %>%
  arrange(pooled_cluster) %>%
  mutate(color = phenotype_palette_master[phenotype_label_short]) %>%
  { stats::setNames(.$color, .$phenotype_label_short) }

p_covariate_clusters <- ggplot(
  cluster_membership_plot_df,
  aes(x = pooled_effect_random, y = y_pos, xmin = pooled_conf_low_random, xmax = pooled_conf_high_random, color = pooled_cluster_label)
) +
  geom_vline(xintercept = 1, linetype = 2, color = "gray45") +
  geom_errorbar(width = 0.10, linewidth = 0.8, orientation = "y") +
  geom_point(size = 2.5) +
  geom_text(aes(x = label_x, label = label_text), hjust = 0, size = 2.8, color = "gray15", show.legend = FALSE) +
  scale_x_log10(
    breaks = c(0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.5, 2.0),
    labels = c("0.6", "0.7", "0.8", "0.9", "1.0", "1.1", "1.25", "1.5", "2.0")
  ) +
  scale_y_continuous(
    breaks = c(4, 3, 2, 1),
    labels = cluster_covariates,
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_color_manual(values = cluster_plot_palette, name = NULL) +
  labs(
    title = "Pooled Covariate Associations With Consensus Cluster Membership",
    x = "Pooled multinomial OR (vs reference cluster, random-effects, log scale)",
    y = NULL
  ) +
  coord_cartesian(xlim = c(0.6, 2.05), clip = "off") +
  theme_manuscript() +
  theme(
    plot.margin = margin(10, 95, 10, 10),
    legend.position = "bottom"
  )

save_plot_safe(p_covariate_clusters, file.path(OUT_DIR, "forest_covariate_effects_cluster_membership.png"), width = 12, height = 10)

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
  filter(site_run %in% model_site_runs) %>%
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

interaction_site_order <- c(sort(unique(cluster_map$site_name[cluster_map$site_run %in% model_site_runs])), "Pooled random-effects")
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
  "Room air" = "#95b46a",
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
  scale_fill_manual(
    values = rs_palette,
    name = "Respiratory support",
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  facet_wrap(~ pooled_cluster_label, ncol = 1) +
  labs(
    title = "Pooled Consortium Respiratory Support Trajectories",
    x = "Hours from ICU admission",
    y = "Patients in cluster-state"
  ) +
  theme_manuscript() +
  theme(
    panel.grid.major.x = element_line(color = "gray90"),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

p_arf_sig <- ggplot(arf_signature_pooled, aes(x = h, y = prop, fill = arf)) +
  geom_area(color = "white", linewidth = 0.15, alpha = 0.98) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = arf_palette,
    name = "ARF subtype",
    guide = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  facet_wrap(~ pooled_cluster_label, ncol = 1) +
  labs(
    title = "Pooled Consortium Acute Respiratory Failure Trajectories",
    x = "Hours from ICU admission",
    y = "Patients in cluster-state"
  ) +
  theme_manuscript() +
  theme(
    panel.grid.major.x = element_line(color = "gray90"),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

save_plot_safe(p_rs_sig, file.path(OUT_DIR, "figure_pooled_resp_support_trajectories.png"), width = 11, height = 14)
save_plot_safe(p_arf_sig, file.path(OUT_DIR, "figure_pooled_arf_trajectories.png"), width = 11, height = 14)

p_trajectory_combo <- cowplot::plot_grid(
  p_rs_sig,
  p_arf_sig,
  labels = c("A", "B"),
  label_size = 16,
  label_fontface = "bold",
  ncol = 2,
  align = "h",
  axis = "tb",
  rel_widths = c(1, 1)
)

save_plot_safe(p_trajectory_combo, file.path(OUT_DIR, "figure_pooled_resp_arf_combo_ab.png"), width = 18, height = 14)

severe_arf_df <- arf_signature_pooled %>%
  mutate(severe_arf = arf %in% c("Hypoxemic ARF", "Mixed ARF")) %>%
  group_by(pooled_cluster, pooled_cluster_label, h) %>%
  summarise(prop = sum(prop[severe_arf], na.rm = TRUE), .groups = "drop")

p_severe_arf <- severe_arf_df %>%
  left_join(distinct(prototype_key, pooled_cluster, phenotype_label_short), by = "pooled_cluster") %>%
  mutate(phenotype_wrapped = stringr::str_wrap(phenotype_label_short, width = 20)) %>%
  ggplot(aes(x = h, y = prop, color = phenotype_wrapped)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(
    values = stats::setNames(unname(cluster_plot_palette_phenotype_short), stringr::str_wrap(names(cluster_plot_palette_phenotype_short), width = 20)),
    name = "Trajectory phenotype"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Severe ARF Burden Over Time by Consensus Trajectory Cluster",
    x = "Hours from ICU admission",
    y = "Share with severe ARF"
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

if (FALSE) {
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

study_years <- 2018:2024

pm25_study_map <- readr::read_csv(file.path(ROOT_DIR, "exposome", "pm25_county_year.csv"), show_col_types = FALSE) %>%
  transmute(
    county_fips = sprintf("%05d", readr::parse_number(GEOID)),
    year = as.integer(year),
    value = as.numeric(pm25_mean)
  ) %>%
  filter(year %in% study_years) %>%
  group_by(county_fips) %>%
  summarise(study_period_mean = mean(value, na.rm = TRUE), .groups = "drop")

no2_study_map <- readr::read_csv(file.path(ROOT_DIR, "exposome", "no2_county_year.csv"), show_col_types = FALSE) %>%
  transmute(
    county_fips = sprintf("%05d", readr::parse_number(GEOID)),
    year = as.integer(year),
    value = as.numeric(no2_mean)
  ) %>%
  filter(year %in% study_years) %>%
  group_by(county_fips) %>%
  summarise(study_period_mean = mean(value, na.rm = TRUE), .groups = "drop")

exposome_map_df <- ggplot2::map_data("county") %>%
  mutate(polyname = paste(region, subregion, sep = ",")) %>%
  left_join(county_fips_lookup, by = "polyname") %>%
  filter(!substr(county_fips, 1, 2) %in% c("02", "15", "72"))

pm25_map_df <- exposome_map_df %>%
  left_join(pm25_study_map, by = "county_fips")

no2_map_df <- exposome_map_df %>%
  left_join(no2_study_map, by = "county_fips")

pm25_limits <- range(pm25_study_map$study_period_mean, na.rm = TRUE)
no2_limits <- range(no2_study_map$study_period_mean, na.rm = TRUE)
pm25_breaks <- scales::breaks_extended(n = 8)(pm25_limits)
no2_breaks <- scales::breaks_extended(n = 8)(no2_limits)

p_pm25_exposome <- ggplot(pm25_map_df, aes(long, lat, group = group, fill = study_period_mean)) +
  geom_polygon(color = NA) +
  coord_quickmap() +
  scale_fill_viridis_c(
    option = "magma",
    limits = pm25_limits,
    breaks = pm25_breaks,
    na.value = "gray94",
    name = expression("Average PM"[2.5] * " (" * mu * "g/m"^3 * "), 2018 to 2024"),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      direction = "horizontal",
      barwidth = unit(10, "cm"),
      barheight = unit(0.5, "cm")
    )
  ) +
  labs(title = "County-Level Average PM₂.₅ Across the Study Period") +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(1.2, "cm"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold")
  )

p_no2_exposome <- ggplot(no2_map_df, aes(long, lat, group = group, fill = study_period_mean)) +
  geom_polygon(color = NA) +
  coord_quickmap() +
  scale_fill_viridis_c(
    option = "plasma",
    limits = no2_limits,
    breaks = no2_breaks,
    na.value = "gray94",
    name = expression("Average NO"[2] * " (ppb), 2018 to 2024"),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      direction = "horizontal",
      barwidth = unit(10, "cm"),
      barheight = unit(0.5, "cm")
    )
  ) +
  labs(title = "County-Level Average NO₂ Across the Study Period") +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(1.2, "cm"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(face = "bold")
  )

write_csv_safe(pm25_study_map, file.path(OUT_DIR, "study_period_pm25_county_mean_2018_2024.csv"))
write_csv_safe(no2_study_map, file.path(OUT_DIR, "study_period_no2_county_mean_2018_2024.csv"))

save_plot_safe(p_pm25_exposome, file.path(OUT_DIR, "map_conus_pm25_study_period_mean.png"), width = 13, height = 8)
save_plot_safe(p_no2_exposome, file.path(OUT_DIR, "map_conus_no2_study_period_mean.png"), width = 13, height = 8)

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
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(0.5, "cm"),
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

p_map_combo <- cowplot::plot_grid(
  p_pm25_exposome,
  p_no2_exposome,
  p_map_dom,
  labels = c("A", "B", "C"),
  label_size = 16,
  label_fontface = "bold",
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1, 1)
)

save_plot_safe(p_map_combo, file.path(OUT_DIR, "figure_maps_exposure_cluster_abc.png"), width = 13, height = 22)
}

consort_flow <- map_dfr(site_payload, function(x) {
  df <- x$consort_flow
  if (is.null(df)) return(tibble())
  df %>%
    mutate(site_run = x$site_run, site_name = x$site_name)
})

consort_reasons <- map_dfr(site_payload, function(x) {
  df <- x$consort_reasons
  if (is.null(df)) return(tibble())
  df %>%
    mutate(site_run = x$site_run, site_name = x$site_name)
})

site_analysis_n <- overall_cont %>%
  filter(variable == "age_years", stratum == "Overall") %>%
  group_by(site_run, site_name) %>%
  summarise(analysis_ready_n = first(n_nonmissing), .groups = "drop")

site_cluster_n <- bycluster_cont %>%
  filter(variable == "age_years", !is.na(stratum)) %>%
  group_by(site_run, site_name) %>%
  summarise(cluster_assigned_n = sum(n_nonmissing, na.rm = TRUE), .groups = "drop")

harmonized_consort_flow <- consort_flow %>%
  filter(step_order <= 4) %>%
  select(site_run, site_name, step_order, step, n_remaining, n_excluded_at_step) %>%
  bind_rows(
    consort_flow %>%
      filter(step_order == 4) %>%
      select(site_run, site_name, prev_n = n_remaining) %>%
      left_join(site_analysis_n, by = c("site_run", "site_name")) %>%
      transmute(
        site_run,
        site_name,
        step_order = 5,
        step = "Lung cancer cases with an ICU stay > 72 hours",
        n_remaining = analysis_ready_n,
        n_excluded_at_step = prev_n - analysis_ready_n
      ),
    site_analysis_n %>%
      left_join(site_cluster_n, by = c("site_run", "site_name")) %>%
      transmute(
        site_run,
        site_name,
        step_order = 6,
        step = "Assigned to respiratory trajectory cluster",
        n_remaining = cluster_assigned_n,
        n_excluded_at_step = analysis_ready_n - cluster_assigned_n
      )
  )

pooled_consort_flow <- harmonized_consort_flow %>%
  group_by(step_order, step) %>%
  summarise(
    n_remaining = sum(n_remaining, na.rm = TRUE),
    n_excluded_at_step = sum(n_excluded_at_step, na.rm = TRUE),
    n_sites = n_distinct(site_name),
    .groups = "drop"
  ) %>%
  arrange(step_order)

pooled_consort_reasons <- consort_reasons %>%
  group_by(reason) %>%
  summarise(
    n = sum(n, na.rm = TRUE),
    n_sites = n_distinct(site_name),
    .groups = "drop"
  ) %>%
  arrange(desc(n))

write_csv_safe(pooled_consort_flow, file.path(OUT_DIR, "publication_table_consort_flow.csv"))
write_csv_safe(pooled_consort_reasons, file.path(OUT_DIR, "publication_table_consort_exclusions.csv"))

wrap_consort_text <- function(x, width = 26) {
  stringr::str_wrap(x, width = width)
}

escape_dot_label <- function(s) {
  s <- as.character(s)
  s <- gsub("\u2212", "-", s, perl = TRUE)
  s <- gsub("\u2013|\u2014", "-", s, perl = TRUE)
  s <- gsub("[\u200B-\u200D\uFEFF]", "", s, perl = TRUE)
  s2 <- suppressWarnings(iconv(s, to = "ASCII//TRANSLIT", sub = ""))
  s2[is.na(s2)] <- s[is.na(s2)]
  s2 <- gsub("\\\\", "\\\\\\\\", s2)
  s2 <- gsub("\"", "\\\\\"", s2)
  s2
}

consort_label_df <- pooled_consort_flow %>%
  mutate(
    step_text = dplyr::recode(
      step,
      "Hospitalizations in date window" = "Hospitalizations in date window",
      ">= 18 years" = "Age 18 years or older",
      "Demographics present" = "Demographics present",
      "Geography present" = "Geography present",
      "Lung cancer dx POA present" = "Lung cancer diagnosis present on admission",
      "Any ICU segment present" = "Any ICU segment present",
      .default = step
    ),
    excl_text = dplyr::case_when(
      step_order == 1 ~ NA_character_,
      step_order == 2 ~ "Excluded: younger than 18 years",
      step_order == 3 ~ "Excluded: missing demographic data",
      step_order == 4 ~ "Excluded: missing county or geographic linkage",
      step_order == 5 ~ "Excluded: no lung cancer ICU stay longer than 72 hours",
      step_order == 6 ~ "Excluded: not assigned to a respiratory trajectory phenotype",
      TRUE ~ "Excluded at this step"
    ),
    main_node = paste0("step", step_order),
    main_label = paste0(wrap_consort_text(step_text, width = 24), "\nN = ", scales::comma(n_remaining)),
    excl_node = paste0("excl", step_order),
    excl_label = ifelse(
      is.na(excl_text),
      NA_character_,
      paste0(wrap_consort_text(excl_text, width = 28), "\nN = ", scales::comma(n_excluded_at_step))
    )
  )

exclusion_df <- consort_label_df %>%
  filter(n_excluded_at_step > 0)

keep_nodes <- vapply(seq_len(nrow(consort_label_df)), function(i) {
  lbl <- sprintf("%s\nN = %s",
                 escape_dot_label(consort_label_df$step_text[i]),
                 scales::comma(consort_label_df$n_remaining[i]))
  sprintf('K%d [label="%s", fillcolor="#F8FBFF", color="#4B5563", fontcolor="#374151", penwidth=1.4];', i, lbl)
}, character(1))

keep_edges <- if (nrow(consort_label_df) > 1) {
  sprintf("K%d -> K%d;", 1:(nrow(consort_label_df) - 1), 2:nrow(consort_label_df))
} else {
  character(0)
}

drop_nodes <- character(0)
drop_edges <- character(0)
rank_same_blocks <- character(0)

for (i in 2:nrow(consort_label_df)) {
  excl_n <- consort_label_df$n_excluded_at_step[i]
  if (is.na(excl_n) || excl_n <= 0) next
  lbl <- sprintf("%s\nN = %s",
                 escape_dot_label(exclusion_df$excl_text[match(i, exclusion_df$step_order)]),
                 scales::comma(excl_n))
  drop_nodes <- c(drop_nodes, sprintf('X%d [label="%s", fillcolor="#FFF7ED", color="#C2410C", fontcolor="#7C2D12", penwidth=1.1];', i, lbl))
  drop_edges <- c(drop_edges, sprintf('K%d -> X%d [color="#C2410C", style=dashed, penwidth=1.0, arrowsize=0.7];', i - 1, i))
  rank_same_blocks <- c(rank_same_blocks, sprintf("{rank=same; K%d; X%d}", i, i))
}

dot <- paste(
  'digraph consort {',
  '  graph [rankdir=TB, nodesep="0.55", ranksep="0.8", margin="0.08", pad="0.2", bgcolor="white"];',
  '  node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=15, margin="0.10,0.08"];',
  '  edge [fontname="Helvetica", color="#6B7280", penwidth=1.2, arrowsize=0.8, arrowhead=normal];',
  paste(keep_nodes, collapse = "\n"),
  paste(drop_nodes, collapse = "\n"),
  paste(keep_edges, collapse = "\n"),
  paste(drop_edges, collapse = "\n"),
  paste(rank_same_blocks, collapse = "\n"),
  '}',
  sep = "\n"
)

consort_diagram <- DiagrammeR::grViz(dot)
save_diagram_safe(consort_diagram, file.path(OUT_DIR, "figure_pooled_consort_diagram.png"), width = 13, height = 9)

cluster_gcs_summary <- map_dfr(site_payload, function(x) {
  df <- x$cluster_gcs_summary
  if (is.null(df)) return(tibble())
  df %>%
    mutate(
      site_run = x$site_run,
      site_name = x$site_name,
      local_cluster = suppressWarnings(as.integer(traj_cluster_ra))
    ) %>%
    filter(!is.na(local_cluster))
})

pooled_sofa_cluster_summary <- static_summary %>%
  left_join(
    cluster_map %>%
      select(site_run, local_cluster, pooled_cluster, pooled_cluster_label, phenotype_label_short, cluster_n = n),
    by = c("site_run", "local_cluster")
  ) %>%
  filter(!is.na(pooled_cluster)) %>%
  group_by(pooled_cluster, pooled_cluster_label, phenotype_label_short) %>%
  summarise(
    n_sites = n_distinct(site_name),
    n_patients = sum(cluster_n, na.rm = TRUE),
    sofa_total_mean = weighted.mean(sofa_total_mean, cluster_n, na.rm = TRUE),
    sofa_resp_mean = weighted.mean(sofa_resp_mean, cluster_n, na.rm = TRUE),
    sofa_cv_mean = weighted.mean(sofa_cv_mean, cluster_n, na.rm = TRUE),
    sofa_renal_mean = weighted.mean(sofa_renal_mean, cluster_n, na.rm = TRUE),
    sofa_liver_mean = weighted.mean(sofa_liver_mean, cluster_n, na.rm = TRUE),
    sofa_coag_mean = weighted.mean(sofa_coag_mean, cluster_n, na.rm = TRUE),
    sofa_cns_mean = weighted.mean(sofa_cns_mean, cluster_n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    cluster_gcs_summary %>%
      left_join(cluster_map %>% select(site_run, local_cluster, pooled_cluster, cluster_n = n), by = c("site_run", "local_cluster")) %>%
      filter(!is.na(pooled_cluster)) %>%
      group_by(pooled_cluster) %>%
      summarise(
        mean_gcs = weighted.mean(mean_gcs, cluster_n, na.rm = TRUE),
        median_gcs_proxy = weighted.mean(median_gcs, cluster_n, na.rm = TRUE),
        pct_gcs_lt8 = weighted.mean(pct_gcs_lt8, cluster_n, na.rm = TRUE),
        pct_gcs_lt13 = weighted.mean(pct_gcs_lt13, cluster_n, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "pooled_cluster"
  )

write_csv_safe(pooled_sofa_cluster_summary, file.path(OUT_DIR, "publication_table_sofa_by_cluster.csv"))

sofa_domain_long <- pooled_sofa_cluster_summary %>%
  select(
    pooled_cluster,
    phenotype_label_short,
    sofa_resp_mean,
    sofa_cv_mean,
    sofa_renal_mean,
    sofa_liver_mean,
    sofa_coag_mean,
    sofa_cns_mean
  ) %>%
  pivot_longer(
    cols = starts_with("sofa_"),
    names_to = "domain",
    values_to = "score"
  ) %>%
  mutate(
    domain = recode(
      domain,
      sofa_resp_mean = "Respiratory",
      sofa_cv_mean = "Cardiovascular",
      sofa_renal_mean = "Renal",
      sofa_liver_mean = "Liver",
      sofa_coag_mean = "Coagulation",
      sofa_cns_mean = "CNS"
    ),
    cluster_label = paste0("C", pooled_cluster, ": ", phenotype_label_short),
    phenotype_label_short = factor(
      phenotype_label_short,
      levels = prototype_key %>%
        distinct(pooled_cluster, phenotype_label_short) %>%
        arrange(pooled_cluster) %>%
        pull(phenotype_label_short)
    ),
    domain = factor(domain, levels = c("Respiratory", "Cardiovascular", "Renal", "Liver", "Coagulation", "CNS"))
  ) %>%
  mutate(
    domain_rev = forcats::fct_rev(domain),
    domain_y = as.numeric(domain_rev),
    phenotype_index = as.integer(phenotype_label_short),
    y = domain_y + c(-0.24, -0.08, 0.08, 0.24)[phenotype_index]
  )

p_sofa_domains <- ggplot(
  sofa_domain_long,
  aes(
    x = score,
    y = y,
    color = phenotype_label_short
  )) +
  geom_segment(
    aes(
      x = 0,
      xend = score,
      y = y,
      yend = y
    ),
    linewidth = 1.2,
    alpha = 0.7
  ) +
  geom_point(
    size = 3.8
  ) +
  geom_text(
    aes(label = sprintf("%.2f", score)),
    nudge_x = 0.035,
    hjust = 0,
    size = 3.2,
    color = "gray20"
  ) +
  scale_color_manual(
    values = cluster_plot_palette_phenotype_short,
    name = "Trajectory phenotype"
  ) +
  scale_x_continuous(
    limits = c(0, max(sofa_domain_long$score, na.rm = TRUE) + 0.45),
    breaks = seq(0, 4, by = 0.5),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    breaks = sort(unique(sofa_domain_long$domain_y)),
    labels = levels(sofa_domain_long$domain_rev),
    expand = expansion(mult = c(0.06, 0.06))
  ) +
  labs(
    title = "Pooled SOFA Domain Signature by Trajectory Phenotype",
    x = "Mean SOFA component score",
    y = NULL
  ) +
  theme_manuscript() +
  theme(
    legend.position = "bottom",
    axis.line.x = element_line(color = "black", linewidth = 0.6),
    axis.ticks.x = element_line(color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 12)
  )

save_plot_safe(p_sofa_domains, file.path(OUT_DIR, "figure_sofa_domain_signature_by_cluster.png"), width = 11.5, height = 6.5)

sofa_totals_long <- pooled_sofa_cluster_summary %>%
  transmute(
    cluster_label = paste0("C", pooled_cluster, ": ", phenotype_label_short),
    pooled_cluster = as.integer(pooled_cluster),
    `Mean total SOFA` = sofa_total_mean,
    `Mean CNS SOFA` = sofa_cns_mean,
    `GCS < 8` = pct_gcs_lt8 * 100,
    `GCS < 13` = pct_gcs_lt13 * 100
  ) %>%
  pivot_longer(
    cols = c(`Mean total SOFA`, `Mean CNS SOFA`, `GCS < 8`, `GCS < 13`),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric, levels = c("Mean total SOFA", "Mean CNS SOFA", "GCS < 8", "GCS < 13")),
    label = case_when(
      metric %in% c("GCS < 8", "GCS < 13") ~ sprintf("%.1f%%", value),
      TRUE ~ sprintf("%.2f", value)
    )
  )

sofa_totals_long <- sofa_totals_long %>%
  mutate(
    metric_group = ifelse(metric %in% c("GCS < 8", "GCS < 13"), "GCS impairment (%)", "SOFA score"),
    cluster_label = factor(cluster_label, levels = paste0("C", prototype_key$pooled_cluster, ": ", prototype_key$phenotype_label_short))
  )

sofa_score_df <- sofa_totals_long %>%
  filter(metric_group == "SOFA score")

gcs_impair_df <- sofa_totals_long %>%
  filter(metric_group == "GCS impairment (%)")

p_sofa_scores <- ggplot(
  sofa_score_df,
  aes(x = metric, y = value, fill = cluster_label)
) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7, color = "white") +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.78),
    vjust = -0.25,
    size = 3.0
  ) +
  scale_fill_manual(values = cluster_plot_palette_short, name = NULL) +
  labs(
    title = "Total SOFA and GCS Impairment by Trajectory Phenotype",
    x = NULL,
    y = "SOFA score"
  ) +
  theme_manuscript() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    axis.ticks.y = element_line(color = "black")
  ) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14)))

p_gcs_impair <- ggplot(
  gcs_impair_df,
  aes(x = metric, y = value, fill = cluster_label)
) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7, color = "white") +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.78),
    vjust = -0.25,
    size = 3.0
  ) +
  scale_fill_manual(values = cluster_plot_palette_short, name = NULL) +
  labs(
    x = NULL,
    y = "Percent with impairment"
  ) +
  theme_manuscript() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    axis.ticks.y = element_line(color = "black")
  ) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14)))

p_sofa_totals <- (p_sofa_scores / p_gcs_impair) +
  plot_layout(heights = c(1, 1), guides = "collect") &
  theme(legend.position = "bottom")

save_plot_safe(p_sofa_totals, file.path(OUT_DIR, "figure_sofa_total_neuro_by_cluster.png"), width = 12, height = 9)

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
    "descriptive_pooling_sites",
    "model_meta_excluded_sites",
    "hopkins_denominator_harmonization",
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
    paste(sort(unique(map_chr(site_payload, "site_name"))), collapse = "; "),
    paste(MODEL_EXCLUDED_SITES, collapse = "; "),
    "Hopkins consort_flow_counts reported 16,958 lung cancer cases with any ICU segment, whereas Hopkins descriptive and cluster exports reflected a larger cohort. For harmonized pooled descriptives, Hopkins count-based descriptive and cluster outputs were rescaled to the consort-flow ICU cohort size while preserving the returned Hopkins means, proportions, and within-cluster distributions.",
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
