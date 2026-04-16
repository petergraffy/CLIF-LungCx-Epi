suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(patchwork)
})

repo_dir <- normalizePath(".", mustWork = TRUE)
pooled_dirs <- list.dirs(file.path(repo_dir, "output"), recursive = FALSE, full.names = TRUE)
pooled_dirs <- pooled_dirs[grepl("pooled_", basename(pooled_dirs))]
stopifnot(length(pooled_dirs) > 0)
latest_pooled_dir <- pooled_dirs[order(basename(pooled_dirs), decreasing = TRUE)][1]

contrast_path <- file.path(latest_pooled_dir, "pooled_modeled_cluster_contrasts.csv")
stopifnot(file.exists(contrast_path))

shorten_phenotype <- function(x) {
  dplyr::recode(
    x,
    "Low-flow oxygen / Minimal ARF" = "Low-flow O2 / Minimal ARF",
    "Low-flow oxygen / Resolving ARF" = "Low-flow O2 / Resolving ARF",
    "Noninvasive ventilation / Resolving ARF" = "NIV / Resolving ARF",
    "Invasive ventilation / Persistent ARF" = "IMV / Persistent ARF",
    .default = x
  )
}

theme_pub <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

save_plot <- function(p, filename, w = 12, h = 9, dpi = 320) {
  ggsave(
    filename = file.path(latest_pooled_dir, filename),
    plot = p,
    width = w,
    height = h,
    dpi = dpi,
    units = "in",
    bg = "white"
  )
}

prototype_path <- file.path(latest_pooled_dir, "pooled_cluster_prototypes.csv")
prototype_tbl <- readr::read_csv(prototype_path, show_col_types = FALSE) %>%
  transmute(
    pooled_cluster = as.integer(pooled_cluster),
    phenotype = phenotype_label
  ) %>%
  arrange(pooled_cluster)

phenotype_palette_master <- c(
  "Minimal ARF" = "#2E8B57",
  "Stable ARF" = "#D4A017",
  "Resolving ARF" = "#2B6CB0",
  "Persistent ARF" = "#C05621"
)
phenotype_palette <- phenotype_palette_master[prototype_tbl$phenotype]

exposure_labels <- c(
  "pm25_5y_z" = "PM2.5",
  "no2_5y_z" = "NO2"
)

contrasts <- readr::read_csv(contrast_path, show_col_types = FALSE) %>%
  left_join(prototype_tbl, by = "pooled_cluster") %>%
  mutate(
    phenotype = factor(phenotype, levels = names(phenotype_palette)),
    pooled_cluster = as.integer(pooled_cluster),
    exposure_label = recode(exposure, !!!exposure_labels)
  ) %>%
  arrange(exposure, pooled_cluster)

dumbbell_df <- contrasts %>%
  transmute(
    exposure_label,
    phenotype,
    low = low_prob_weighted,
    high = high_prob_weighted,
    risk_diff_weighted,
    risk_ratio_weighted,
    n_sites,
    n_patients
  )

p_dumbbell <- ggplot(dumbbell_df, aes(y = phenotype)) +
  geom_segment(aes(x = low, xend = high, yend = phenotype), linewidth = 2.2, color = "grey55", lineend = "round") +
  geom_point(aes(x = low), size = 3.2, shape = 21, stroke = 0.7, fill = "white", color = "grey25") +
  geom_point(aes(x = high, fill = phenotype), size = 3.4, shape = 21, stroke = 0.7, color = "white") +
  geom_text(
    aes(
      x = pmax(low, high) + 0.012,
      label = sprintf("%+.1f%%", 100 * risk_diff_weighted)
    ),
    hjust = 0,
    size = 3.3,
    color = "grey20"
  ) +
  scale_fill_manual(values = phenotype_palette, drop = FALSE) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0.04, 0.56),
    breaks = seq(0.05, 0.55, by = 0.05),
    expand = expansion(mult = c(0.01, 0.14))
  ) +
  facet_wrap(~ exposure_label, ncol = 1) +
  labs(
    title = "Pooled effect modification mockup: modeled mortality at low vs high exposure",
    x = "Predicted probability of death or hospice discharge",
    y = NULL,
    fill = "Phenotype"
  ) +
  annotate("text", x = 0.045, y = length(levels(dumbbell_df$phenotype)) + 0.55, label = "Open circle = low exposure\nFilled circle = high exposure", hjust = 0, vjust = 1, size = 3.2) +
  theme_pub() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12),
    plot.margin = margin(10, 20, 10, 10)
  )

risk_ratio_df <- contrasts %>%
  transmute(
    exposure_label,
    phenotype,
    risk_ratio_weighted,
    risk_diff_weighted
  )

p_rr <- ggplot(risk_ratio_df, aes(x = risk_ratio_weighted, y = phenotype, color = phenotype)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey70") +
  geom_point(size = 3.1) +
  geom_text(
    aes(
      x = risk_ratio_weighted,
      label = sprintf("RR %.02f", risk_ratio_weighted)
    ),
    nudge_x = 0.012,
    hjust = 0,
    size = 3.2,
    color = "grey20"
  ) +
  scale_color_manual(values = phenotype_palette, drop = FALSE) +
  scale_x_continuous(
    limits = c(0.88, 1.16),
    breaks = seq(0.90, 1.15, by = 0.05),
    expand = expansion(mult = c(0.01, 0.18))
  ) +
  facet_wrap(~ exposure_label, ncol = 1) +
  labs(
    title = "Relative mortality contrast at high vs low exposure",
    x = "Risk ratio",
    y = NULL,
    color = "Phenotype"
  ) +
  theme_pub() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12)
  )

p_combo <- p_dumbbell / p_rr + plot_layout(heights = c(2.2, 1.2))

save_plot(p_combo, "figure_effect_modification_mockup.png", w = 12, h = 11)
save_plot(p_dumbbell, "figure_effect_modification_mockup_dumbbell.png", w = 11, h = 8.5)
save_plot(p_rr, "figure_effect_modification_mockup_rr.png", w = 11, h = 7)

readr::write_csv(dumbbell_df, file.path(latest_pooled_dir, "effect_modification_mockup_table.csv"))

cat("Saved effect modification mockups to: ", latest_pooled_dir, "\n", sep = "")
