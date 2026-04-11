suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(patchwork)
})

repo_dir <- normalizePath(".", mustWork = TRUE)
pooled_dirs <- list.dirs(file.path(repo_dir, "output"), recursive = FALSE, full.names = TRUE)
pooled_dirs <- pooled_dirs[grepl("pooled_", basename(pooled_dirs))]
stopifnot(length(pooled_dirs) > 0)
latest_pooled_dir <- pooled_dirs[order(file.info(pooled_dirs)$mtime, decreasing = TRUE)][1]

contrast_path <- file.path(latest_pooled_dir, "pooled_modeled_cluster_contrasts.csv")
stopifnot(file.exists(contrast_path))

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

phenotype_palette <- c(
  "Room air / No ARF" = "#2E8B57",
  "Low-flow O2 / No ARF" = "#D4A017",
  "NIV / No ARF" = "#2B6CB0",
  "IMV / Predom. no ARF" = "#C05621",
  "IMV / Hypoxemic/mixed ARF" = "#B83280"
)

exposure_labels <- c(
  "pm25_5y_z" = "PM2.5",
  "no2_5y_z" = "NO2"
)

contrasts <- readr::read_csv(contrast_path, show_col_types = FALSE) %>%
  mutate(
    phenotype = stringr::str_remove(pooled_cluster_label, "^Cluster \\d+: "),
    phenotype = recode(
      phenotype,
      "Low-flow oxygen / No ARF" = "Low-flow O2 / No ARF",
      "Noninvasive ventilation / No ARF" = "NIV / No ARF",
      "Invasive ventilation / Predominantly no ARF" = "IMV / Predom. no ARF",
      "Invasive ventilation / Hypoxemic or mixed ARF" = "IMV / Hypoxemic/mixed ARF"
    ),
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
  geom_segment(aes(x = low, xend = high, yend = phenotype), linewidth = 1.5, color = "grey70") +
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
    limits = c(0.04, 0.48),
    breaks = seq(0.05, 0.45, by = 0.05),
    expand = expansion(mult = c(0.01, 0.14))
  ) +
  facet_wrap(~ exposure_label, ncol = 1) +
  labs(
    title = "Pooled effect modification mockup: modeled mortality at low vs high exposure",
    x = "Predicted probability of death or hospice discharge",
    y = NULL,
    fill = "Phenotype"
  ) +
  annotate("text", x = 0.045, y = 5.55, label = "Open circle = low exposure\nFilled circle = high exposure", hjust = 0, vjust = 1, size = 3.2) +
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
