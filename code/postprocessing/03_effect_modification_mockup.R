suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(patchwork)
  library(cowplot)
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
  "pm25_5y_z" = "PM\u2082.\u2085",
  "no2_5y_z" = "NO\u2082"
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
  geom_point(aes(x = low), size = 5.8, shape = 21, stroke = 1.0, fill = "white", color = "grey25") +
  geom_point(aes(x = high, fill = phenotype), size = 6.0, shape = 21, stroke = 1.0, color = "white") +
  geom_text(
    aes(
      x = pmax(low, high) + 0.012,
      label = sprintf("%+.1f%%", 100 * risk_diff_weighted)
    ),
    hjust = 0,
    size = 4.1,
    color = "grey20"
  ) +
  scale_fill_manual(values = phenotype_palette, drop = FALSE) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0.04, 0.56),
    breaks = seq(0.05, 0.55, by = 0.05),
    expand = expansion(mult = c(0.01, 0.14))
  ) +
  facet_wrap(~ exposure_label, ncol = 1, scales = "fixed") +
  labs(
    title = "Adjusted death or hospice probability\nat low versus high air pollution exposure by phenotype",
    subtitle = "Low exposure = 1 SD below the site mean; high exposure = 1 SD above the site mean.\nOpen circles indicate low exposure and filled circles indicate high exposure.",
    x = "Predicted probability of death or hospice discharge",
    y = NULL,
    fill = "Phenotype"
  ) +
  theme_pub(base_size = 15) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    plot.margin = margin(10, 12, 10, 10)
  )

risk_diff_df <- contrasts %>%
  transmute(
    exposure_label,
    phenotype,
    risk_diff_weighted = 100 * (high_prob_weighted - low_prob_weighted)
  )

p_rr <- ggplot(risk_diff_df, aes(x = risk_diff_weighted, y = phenotype, color = phenotype)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey70") +
  geom_segment(aes(x = 0, xend = risk_diff_weighted, yend = phenotype), linewidth = 1.8, alpha = 0.45, lineend = "round") +
  geom_point(size = 4.2) +
  geom_text(
    aes(
      x = risk_diff_weighted,
      label = sprintf("%+.1f%%", risk_diff_weighted)
    ),
    nudge_x = 0.012,
    hjust = 0,
    size = 3.2,
    color = "grey20"
  ) +
  scale_color_manual(values = phenotype_palette, drop = FALSE) +
  scale_x_continuous(
    limits = c(-8, 10),
    breaks = c(-6, -3, 0, 3, 6, 9),
    labels = function(x) sprintf("%+.0f%%", x),
    expand = expansion(mult = c(0.01, 0.12))
  ) +
  facet_wrap(~ exposure_label, ncol = 1) +
  labs(
    x = "Absolute risk difference at high vs low exposure",
    y = NULL,
    color = "Phenotype"
  ) +
  theme_pub() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12),
    plot.margin = margin(8, 8, 8, 0)
  )

p_title <- ggdraw() +
  draw_label(
    "Pooled mortality contrasts at low versus high long-term air pollution exposure across trajectory phenotypes",
    fontface = "bold",
    x = 0,
    hjust = 0,
    size = 14
  )

p_combo_body <- p_dumbbell + p_rr + plot_layout(widths = c(1.65, 1), guides = "collect")
p_combo <- (p_title / p_combo_body + plot_layout(heights = c(0.08, 1))) &
  theme(legend.position = "bottom")

save_plot(p_combo, "figure_effect_modification_mockup.png", w = 15, h = 8.5)
save_plot(p_dumbbell, "figure_effect_modification_mockup_dumbbell.png", w = 13.5, h = 9.5)
save_plot(p_rr, "figure_effect_modification_mockup_rr.png", w = 8.5, h = 8.5)

readr::write_csv(dumbbell_df, file.path(latest_pooled_dir, "effect_modification_mockup_table.csv"))

cat("Saved effect modification mockups to: ", latest_pooled_dir, "\n", sep = "")
