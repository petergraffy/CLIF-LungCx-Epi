suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

ROOT_DIR <- "."
SITE_DIR <- file.path(ROOT_DIR, "sites")

pooled_dirs <- list.dirs(file.path(ROOT_DIR, "output"), recursive = FALSE, full.names = TRUE)
pooled_dirs <- pooled_dirs[grepl("pooled_[0-9]{8}$", basename(pooled_dirs))]
stopifnot(length(pooled_dirs) > 0)
OUT_DIR <- pooled_dirs[which.max(as.integer(sub("pooled_", "", basename(pooled_dirs))))]

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

fmt_effect <- function(est, low, high, digits = 2) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"), est, low, high)
}

meta_random_effects <- function(df) {
  yi <- df$estimate
  sei <- df$std_error
  keep <- is.finite(yi) & is.finite(sei) & sei > 0
  yi <- yi[keep]
  sei <- sei[keep]
  k <- length(yi)
  if (k == 0) {
    return(tibble(
      k = 0,
      estimate_random = NA_real_,
      std_error_random = NA_real_,
      conf_low_random = NA_real_,
      conf_high_random = NA_real_,
      tau2 = NA_real_,
      q = NA_real_,
      i2 = NA_real_
    ))
  }

  wi <- 1 / (sei ^ 2)
  mu_fixed <- sum(wi * yi) / sum(wi)
  q <- sum(wi * (yi - mu_fixed) ^ 2)
  c_val <- sum(wi) - (sum(wi ^ 2) / sum(wi))
  tau2 <- if (k > 1 && c_val > 0) max((q - (k - 1)) / c_val, 0) else 0
  wi_re <- 1 / (sei ^ 2 + tau2)
  mu_re <- sum(wi_re * yi) / sum(wi_re)
  se_re <- sqrt(1 / sum(wi_re))

  tibble(
    k = k,
    estimate_random = mu_re,
    std_error_random = se_re,
    conf_low_random = mu_re - 1.96 * se_re,
    conf_high_random = mu_re + 1.96 * se_re,
    tau2 = tau2,
    q = q,
    i2 = if (k > 1 && q > 0) max((q - (k - 1)) / q, 0) else 0
  )
}

phenotype_palette <- c(
  "Resolving ARF" = "#2B6CB0",
  "Stable ARF" = "#D4A017",
  "Persistent ARF" = "#C05621"
)

exposure_palette <- c(
  "PM2.5" = "#4C78A8",
  "NO2" = "#E4572E"
)

site_files <- list.files(SITE_DIR, pattern = "model_outputs_standardized_.*\\.csv$", recursive = TRUE, full.names = TRUE)
stopifnot(length(site_files) > 0)

cluster_map <- readr::read_csv(file.path(OUT_DIR, "pooled_cluster_mapping.csv"), show_col_types = FALSE) %>%
  transmute(
    site_run,
    local_cluster = as.integer(local_cluster),
    pooled_cluster = as.integer(pooled_cluster),
    pooled_cluster_label = as.character(pooled_cluster_label)
  )

prototype_key <- readr::read_csv(file.path(OUT_DIR, "pooled_cluster_prototypes.csv"), show_col_types = FALSE) %>%
  transmute(
    pooled_cluster = as.integer(pooled_cluster),
    pooled_cluster_label = as.character(pooled_cluster_label),
    phenotype_label_short = str_replace(pooled_cluster_label, "^Cluster [0-9]+: ", "")
  )

multinom_terms <- purrr::map_dfr(site_files, function(path) {
  site_run <- basename(dirname(path))
  site_name <- site_run
  readr::read_csv(path, show_col_types = FALSE) %>%
    filter(model == "multinom_traj_cluster", term %in% c("pm25_5y_z", "no2_5y_z")) %>%
    mutate(
      site_run = site_run,
      site_name = site_name,
      local_cluster = suppressWarnings(as.integer(outcome_level))
    ) %>%
    left_join(cluster_map, by = c("site_run", "local_cluster")) %>%
    filter(!is.na(pooled_cluster), pooled_cluster != 1L) %>%
    transmute(
      site_run,
      site_name,
      pooled_cluster,
      pooled_cluster_label,
      term,
      estimate,
      std_error,
      effect,
      conf_low,
      conf_high,
      p_value
    )
})

pooled_cluster_membership_pollution <- multinom_terms %>%
  group_by(pooled_cluster, pooled_cluster_label, term) %>%
  group_modify(~ {
    meta_random_effects(.x) %>%
      mutate(site_included = paste(sort(unique(.x$site_name)), collapse = "; "))
  }) %>%
  ungroup() %>%
  left_join(prototype_key, by = c("pooled_cluster", "pooled_cluster_label")) %>%
  mutate(
    exposure = recode(term, pm25_5y_z = "PM2.5", no2_5y_z = "NO2"),
    pooled_effect_random = exp(estimate_random),
    pooled_conf_low_random = exp(conf_low_random),
    pooled_conf_high_random = exp(conf_high_random),
    estimate_ci = fmt_effect(pooled_effect_random, pooled_conf_low_random, pooled_conf_high_random),
    phenotype_label_short = factor(
      phenotype_label_short,
      levels = c("Resolving ARF", "Stable ARF", "Persistent ARF")
    ),
    exposure = factor(exposure, levels = c("NO2", "PM2.5"))
  ) %>%
  arrange(exposure, phenotype_label_short)

reference_rows <- crossing(
  pooled_cluster = 1L,
  pooled_cluster_label = "Cluster 1: Minimal ARF",
  term = c("no2_5y_z", "pm25_5y_z"),
  exposure = c("NO2", "PM2.5")
) %>%
  mutate(
    k = 10L,
    estimate_random = 0,
    std_error_random = NA_real_,
    conf_low_random = 0,
    conf_high_random = 0,
    tau2 = NA_real_,
    q = NA_real_,
    i2 = NA_real_,
    site_included = paste(sort(unique(multinom_terms$site_name)), collapse = "; "),
    phenotype_label_short = factor("Minimal ARF", levels = c("Minimal ARF", "Resolving ARF", "Stable ARF", "Persistent ARF")),
    pooled_effect_random = 1,
    pooled_conf_low_random = 1,
    pooled_conf_high_random = 1,
    estimate_ci = "Reference"
  )

pooled_cluster_membership_pollution_plot <- bind_rows(
  pooled_cluster_membership_pollution,
  reference_rows
) %>%
  mutate(
    phenotype_label_short = factor(
      as.character(phenotype_label_short),
      levels = c("Minimal ARF", "Resolving ARF", "Stable ARF", "Persistent ARF")
    ),
    exposure = factor(exposure, levels = c("NO2", "PM2.5"))
  )

write_csv_safe(
  pooled_cluster_membership_pollution,
  file.path(OUT_DIR, "pooled_cluster_membership_pollution.csv")
)

plot_df <- pooled_cluster_membership_pollution_plot %>%
  mutate(
    phenotype_label_short = forcats::fct_rev(phenotype_label_short),
    y_base = c(
      "Persistent ARF" = 4,
      "Stable ARF" = 3,
      "Resolving ARF" = 2,
      "Minimal ARF" = 1
    )[as.character(phenotype_label_short)],
    y_pos = y_base + c("NO2" = 0.12, "PM2.5" = -0.12)[as.character(exposure)],
    label_x = case_when(
      estimate_ci == "Reference" ~ 1.28,
      TRUE ~ pmin(pooled_conf_high_random * 1.05, 2.15)
    )
  )

p_cluster_membership_pollution <- ggplot(
  plot_df,
  aes(
    x = pooled_effect_random,
    y = y_pos,
    xmin = pooled_conf_low_random,
    xmax = pooled_conf_high_random,
    color = exposure
  )
) +
  geom_vline(xintercept = 1, linetype = 2, color = "gray45") +
  geom_errorbar(width = 0.12, linewidth = 0.9, orientation = "y") +
  geom_point(size = 3) +
  geom_text(
    aes(x = label_x, label = estimate_ci),
    hjust = 0,
    size = 3.2,
    color = "gray15",
    show.legend = FALSE
  ) +
  scale_x_log10(
    breaks = c(0.6, 0.8, 1.0, 1.2, 1.5, 2.0),
    labels = c("0.6", "0.8", "1.0", "1.2", "1.5", "2.0")
  ) +
  scale_y_continuous(
    breaks = c(1, 2, 3, 4),
    labels = c("Minimal ARF", "Resolving ARF", "Stable ARF", "Persistent ARF"),
    expand = expansion(mult = c(0.08, 0.08))
  ) +
  scale_color_manual(
    values = exposure_palette,
    breaks = c("NO2", "PM2.5"),
    labels = c(expression(NO[2]), expression(PM[2.5])),
    name = NULL
  ) +
  labs(
    title = "Pooled Air Pollution Associations With Respiratory Phenotype Membership",
    x = "Pooled multinomial OR (log scale)",
    y = NULL
  ) +
  coord_cartesian(xlim = c(0.62, 2.2), clip = "off") +
  theme_manuscript() +
  theme(
    plot.margin = margin(12, 120, 12, 12),
    legend.position = "bottom",
    plot.title = element_text(size = 18, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 12),
    legend.text = element_text(size = 13),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.8)
  )

save_plot_safe(
  p_cluster_membership_pollution,
  file.path(OUT_DIR, "forest_air_pollution_cluster_membership.png"),
  width = 12.5,
  height = 6.5
)

message("Cluster-membership pollution forest written to: ", OUT_DIR)
