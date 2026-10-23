suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
  library(patchwork)
})

repo_dir <- normalizePath(".", mustWork = TRUE)
pooled_dirs <- list.dirs(file.path(repo_dir, "output"), recursive = FALSE, full.names = TRUE)
pooled_dirs <- pooled_dirs[grepl("pooled_[0-9]{8}$", basename(pooled_dirs))]
stopifnot(length(pooled_dirs) > 0)
latest_pooled_dir <- pooled_dirs[order(basename(pooled_dirs), decreasing = TRUE)][1]

save_plot <- function(p, filename, w = 12, h = 8, dpi = 320) {
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

theme_pub <- function(base_size = 15) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.line.x = element_line(color = "black", linewidth = 0.6),
      axis.line.y = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.55)
    )
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

phenotype_levels <- c("Resolved", "Mild", "Moderate", "Persistent")
phenotype_palette <- c(
  "Resolved" = "#2E8B57",
  "Mild" = "#2B6CB0",
  "Moderate" = "#D4A017",
  "Persistent" = "#C05621"
)

prototype_tbl <- readr::read_csv(file.path(latest_pooled_dir, "pooled_cluster_prototypes.csv"), show_col_types = FALSE) %>%
  transmute(
    pooled_cluster = as.integer(pooled_cluster),
    phenotype = factor(phenotype_label_short, levels = phenotype_levels),
    cluster_label = paste0("C", pooled_cluster, ": ", phenotype_label_short)
  ) %>%
  arrange(pooled_cluster)

by_cluster_cat <- readr::read_csv(file.path(latest_pooled_dir, "pooled_table1_by_cluster_categorical.csv"), show_col_types = FALSE) %>%
  mutate(pooled_cluster = as.integer(pooled_cluster)) %>%
  left_join(prototype_tbl, by = "pooled_cluster")

by_cluster_cont <- readr::read_csv(file.path(latest_pooled_dir, "pooled_table1_by_cluster_continuous.csv"), show_col_types = FALSE) %>%
  mutate(pooled_cluster = as.integer(pooled_cluster)) %>%
  left_join(prototype_tbl, by = "pooled_cluster")

cluster_map <- readr::read_csv(file.path(latest_pooled_dir, "pooled_cluster_mapping.csv"), show_col_types = FALSE) %>%
  transmute(
    site_run,
    local_cluster = as.integer(local_cluster),
    pooled_cluster = as.integer(pooled_cluster),
    phenotype = factor(phenotype_label_short, levels = phenotype_levels)
  )

model_output_files <- list.files(
  file.path(repo_dir, "sites"),
  pattern = "^model_outputs_standardized.*\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

model_outputs <- purrr::map_dfr(model_output_files, function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(site_run = basename(dirname(path)))
})

binary_outcome_labels <- c(
  "death_or_hospice" = "Death or hospice",
  "death_in_hosp" = "In-hospital death",
  "hospice_discharge" = "Hospice discharge",
  "any_imv_72h" = "Any IMV in first 72h",
  "vaso_any_72h" = "Any vasopressor in first 72h"
)

binary_outcomes <- by_cluster_cat %>%
  filter(variable %in% names(binary_outcome_labels), level == "1") %>%
  transmute(
    pooled_cluster,
    phenotype,
    outcome = recode(variable, !!!binary_outcome_labels),
    events = as.numeric(n),
    total = as.numeric(denom),
    risk = as.numeric(prop)
  ) %>%
  mutate(
    phenotype = factor(phenotype, levels = phenotype_levels),
    outcome = factor(outcome, levels = unname(binary_outcome_labels))
  )

continuous_outcome_labels <- c(
  "sofa_total" = "SOFA total",
  "icu_los_hours" = "ICU length of stay, days",
  "imv_hours_72h" = "IMV hours in first 72h",
  "vaso_hours_72h" = "Vasopressor hours in first 72h"
)

continuous_outcomes <- by_cluster_cont %>%
  filter(variable %in% names(continuous_outcome_labels)) %>%
  transmute(
    pooled_cluster,
    phenotype,
    outcome = recode(variable, !!!continuous_outcome_labels),
    n = as.numeric(n_nonmissing),
    mean = if_else(variable == "icu_los_hours", as.numeric(mean) / 24, as.numeric(mean)),
    sd = if_else(variable == "icu_los_hours", as.numeric(sd) / 24, as.numeric(sd)),
    variable
  ) %>%
  mutate(
    phenotype = factor(phenotype, levels = phenotype_levels),
    outcome = factor(outcome, levels = unname(continuous_outcome_labels))
  )

readr::write_csv(binary_outcomes, file.path(latest_pooled_dir, "publication_table_phenotype_binary_outcomes.csv"))
readr::write_csv(continuous_outcomes, file.path(latest_pooled_dir, "publication_table_phenotype_continuous_outcomes.csv"))

binary_matrix <- binary_outcomes %>%
  mutate(
    display_value = percent(risk, accuracy = 1),
    row_group = "Clinical outcomes",
    row_label = as.character(outcome),
    fill_value = risk
  ) %>%
  select(row_group, row_label, phenotype, display_value, fill_value)

continuous_matrix <- continuous_outcomes %>%
  mutate(
    display_value = number(mean, accuracy = 0.1),
    row_group = "Organ support",
    row_label = as.character(outcome),
    fill_value = scales::rescale(mean, to = c(0, 1))
  ) %>%
  select(row_group, row_label, phenotype, display_value, fill_value)

outcome_matrix <- bind_rows(binary_matrix, continuous_matrix) %>%
  mutate(
    row_label = factor(
      row_label,
      levels = rev(c(
        "Death or hospice",
        "In-hospital death",
        "Hospice discharge",
        "Any IMV in first 72h",
        "Any vasopressor in first 72h",
        "SOFA total",
        "ICU length of stay, days",
        "IMV hours in first 72h",
        "Vasopressor hours in first 72h"
      ))
    ),
    phenotype = factor(phenotype, levels = phenotype_levels)
  )

p_outcome_matrix <- ggplot(outcome_matrix, aes(x = phenotype, y = row_label, fill = fill_value)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = display_value), size = 4.7, fontface = "bold", color = "gray12") +
  facet_grid(row_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient(low = "#F4F7EC", high = "#C05621", guide = "none") +
  labs(
    title = "Phenotype Outcome and Organ Support Burden",
    x = NULL,
    y = NULL
  ) +
  theme_pub(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 22, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold", angle = 20, hjust = 1),
    axis.text.y = element_text(size = 14),
    strip.text.y.left = element_text(size = 14, face = "bold", angle = 90),
    strip.background = element_rect(fill = "#F4F6F8", color = "gray70"),
    strip.placement = "outside"
  )

save_plot(p_outcome_matrix, "figure_phenotype_outcome_matrix.png", w = 12.5, h = 7.8)

p_burden_binary <- binary_outcomes %>%
  mutate(outcome = fct_rev(outcome)) %>%
  ggplot(aes(x = risk, y = outcome, color = phenotype)) +
  geom_line(aes(group = outcome), color = "gray82", linewidth = 1.3) +
  geom_point(size = 4.5) +
  geom_text(
    aes(label = percent(risk, accuracy = 1)),
    nudge_x = 0.025,
    size = 4.1,
    color = "gray18",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = phenotype_palette,
    breaks = c("Mild", "Moderate", "Persistent")
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.01, 0.08))
  ) +
  labs(
    title = "Clinical Outcomes by Respiratory Trajectory Phenotype",
    x = "Patients with outcome",
    y = NULL,
    color = "Phenotype"
  ) +
  theme_pub(base_size = 16) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 17, face = "bold"),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15, face = "bold")
  )

p_burden_continuous <- continuous_outcomes %>%
  mutate(outcome = fct_rev(outcome)) %>%
  ggplot(aes(x = mean, y = outcome, color = phenotype)) +
  geom_line(aes(group = outcome), color = "gray82", linewidth = 1.3) +
  geom_point(size = 4.5) +
  geom_text(
    aes(label = number(mean, accuracy = 0.1)),
    nudge_x = 0.7,
    size = 4.1,
    color = "gray18",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = phenotype_palette,
    breaks = c("Mild", "Moderate", "Persistent")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.15))) +
  labs(
    title = "Organ Dysfunction and Support Burden by Phenotype",
    x = "Mean value",
    y = NULL,
    color = "Phenotype"
  ) +
  theme_pub(base_size = 16) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 17, face = "bold"),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15, face = "bold")
  )

p_outcome_burden <- p_burden_binary / p_burden_continuous +
  plot_layout(guides = "collect", heights = c(1.05, 1)) &
  theme(legend.position = "bottom")

save_plot(p_outcome_burden, "figure_phenotype_outcome_burden.png", w = 13.5, h = 11)
save_plot(p_burden_binary, "figure_phenotype_binary_outcomes.png", w = 12.5, h = 6.8)
save_plot(p_burden_continuous, "figure_phenotype_continuous_burden.png", w = 12.5, h = 6.8)

ref_binary <- binary_outcomes %>%
  filter(phenotype == "Resolved") %>%
  select(outcome, ref_events = events, ref_total = total, ref_risk = risk)

binary_rr <- binary_outcomes %>%
  filter(phenotype != "Resolved") %>%
  left_join(ref_binary, by = "outcome") %>%
  mutate(
    rr = risk / ref_risk,
    se_log_rr = sqrt((1 / events) - (1 / total) + (1 / ref_events) - (1 / ref_total)),
    conf_low = exp(log(rr) - 1.96 * se_log_rr),
    conf_high = exp(log(rr) + 1.96 * se_log_rr),
    estimate_ci = sprintf("%.2f (%.2f, %.2f)", rr, conf_low, conf_high)
  )

readr::write_csv(binary_rr, file.path(latest_pooled_dir, "publication_table_phenotype_binary_rr_vs_minimal.csv"))

p_binary_rr <- binary_rr %>%
  mutate(outcome = fct_rev(outcome)) %>%
  ggplot(aes(x = rr, y = outcome, xmin = conf_low, xmax = conf_high, color = phenotype)) +
  geom_vline(xintercept = 1, linetype = 2, color = "gray45") +
  geom_errorbar(width = 0.18, linewidth = 0.9, orientation = "y", position = position_dodge(width = 0.55)) +
  geom_point(size = 3.9, position = position_dodge(width = 0.55)) +
  scale_color_manual(
    values = phenotype_palette,
    breaks = c("Mild", "Moderate", "Persistent")
  ) +
  scale_x_log10(
    breaks = c(0.75, 1, 1.5, 2, 3, 5, 8, 12, 16, 20),
    labels = c("0.75", "1.0", "1.5", "2", "3", "5", "8", "12", "16", "20")
  ) +
  coord_cartesian(xlim = c(0.75, 20), clip = "off") +
  labs(
    title = "Outcome Risk by Phenotype Compared With Resolved",
    x = "Risk ratio vs Resolved (log scale)",
    y = NULL,
    color = "Phenotype"
  ) +
  theme_pub(base_size = 16) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 21, face = "bold"),
    axis.title = element_text(size = 17, face = "bold"),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15, face = "bold"),
    plot.margin = margin(10, 20, 10, 10)
  )

save_plot(p_binary_rr, "forest_phenotype_binary_outcomes_vs_minimal.png", w = 13.5, h = 7)

ref_cont <- continuous_outcomes %>%
  filter(phenotype == "Resolved") %>%
  select(outcome, ref_n = n, ref_mean = mean, ref_sd = sd)

continuous_md <- continuous_outcomes %>%
  filter(phenotype != "Resolved") %>%
  left_join(ref_cont, by = "outcome") %>%
  mutate(
    mean_diff = mean - ref_mean,
    se_diff = sqrt((sd^2 / n) + (ref_sd^2 / ref_n)),
    conf_low = mean_diff - 1.96 * se_diff,
    conf_high = mean_diff + 1.96 * se_diff,
    estimate_ci = sprintf("%.1f (%.1f, %.1f)", mean_diff, conf_low, conf_high)
  )

readr::write_csv(continuous_md, file.path(latest_pooled_dir, "publication_table_phenotype_continuous_mean_difference_vs_minimal.csv"))

p_continuous_md <- continuous_md %>%
  mutate(
    phenotype = factor(phenotype, levels = rev(c("Mild", "Moderate", "Persistent"))),
    outcome = factor(outcome, levels = unname(continuous_outcome_labels))
  ) %>%
  ggplot(aes(x = mean_diff, y = phenotype, xmin = conf_low, xmax = conf_high, color = phenotype)) +
  geom_vline(xintercept = 0, linetype = 2, color = "gray45") +
  geom_errorbar(width = 0.18, linewidth = 0.9, orientation = "y") +
  geom_point(size = 4.2) +
  facet_wrap(~ outcome, scales = "free_x", ncol = 2) +
  scale_color_manual(values = phenotype_palette, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(
    title = "Mean Difference in Organ Support Burden Compared With Resolved",
    x = "Mean difference vs Resolved",
    y = NULL,
    color = "Phenotype"
  ) +
  theme_pub(base_size = 16) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 21, face = "bold"),
    axis.title = element_text(size = 17, face = "bold"),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15, face = "bold")
  )

save_plot(p_continuous_md, "forest_phenotype_continuous_outcomes_vs_minimal.png", w = 13.5, h = 8)

phenotype_model_terms <- model_outputs %>%
  filter(
    model %in% c("glm_mortality", "glm_los_quasipoisson", "cox72_landmark"),
    str_detect(term, "^traj_cluster_ra[0-9]+$")
  ) %>%
  mutate(
    local_cluster = as.integer(str_remove(term, "^traj_cluster_ra")),
    model_outcome = recode(
      model,
      "glm_mortality" = "Death or hospice",
      "glm_los_quasipoisson" = "ICU length of stay",
      "cox72_landmark" = "Post-72h death or hospice"
    ),
    estimand = recode(
      model,
      "glm_mortality" = "Adjusted odds ratio",
      "glm_los_quasipoisson" = "Adjusted mean ratio",
      "cox72_landmark" = "Adjusted hazard ratio"
    )
  ) %>%
  left_join(cluster_map, by = c("site_run", "local_cluster")) %>%
  filter(!is.na(pooled_cluster), pooled_cluster != 1L)

phenotype_model_effects <- phenotype_model_terms %>%
  group_by(model_outcome, estimand, model, effect_type, pooled_cluster, phenotype) %>%
  group_modify(~ {
    meta_random_effects(.x) %>%
      mutate(site_included = paste(sort(unique(.x$site_run)), collapse = "; "))
  }) %>%
  ungroup() %>%
  mutate(
    pooled_effect_random = exp(estimate_random),
    pooled_conf_low_random = exp(conf_low_random),
    pooled_conf_high_random = exp(conf_high_random),
    estimate_ci = sprintf("%.2f (%.2f, %.2f)", pooled_effect_random, pooled_conf_low_random, pooled_conf_high_random)
  ) %>%
  arrange(model_outcome, pooled_cluster)

phenotype_model_reference <- tidyr::expand_grid(
  model_outcome = c("Death or hospice", "Post-72h death or hospice", "ICU length of stay"),
  phenotype = factor("Resolved", levels = phenotype_levels)
) %>%
  mutate(
    estimand = case_when(
      model_outcome == "Death or hospice" ~ "Adjusted odds ratio",
      model_outcome == "Post-72h death or hospice" ~ "Adjusted hazard ratio",
      TRUE ~ "Adjusted mean ratio"
    ),
    model = case_when(
      model_outcome == "Death or hospice" ~ "glm_mortality",
      model_outcome == "Post-72h death or hospice" ~ "cox72_landmark",
      TRUE ~ "glm_los_quasipoisson"
    ),
    effect_type = case_when(
      model_outcome == "Death or hospice" ~ "OR",
      model_outcome == "Post-72h death or hospice" ~ "HR",
      TRUE ~ "RR_like"
    ),
    pooled_cluster = 1L,
    k = n_distinct(cluster_map$site_run),
    estimate_random = 0,
    std_error_random = NA_real_,
    conf_low_random = 0,
    conf_high_random = 0,
    tau2 = NA_real_,
    q = NA_real_,
    i2 = NA_real_,
    site_included = paste(sort(unique(cluster_map$site_run)), collapse = "; "),
    pooled_effect_random = 1,
    pooled_conf_low_random = 1,
    pooled_conf_high_random = 1,
    estimate_ci = "Reference"
  )

phenotype_model_plot_df <- bind_rows(phenotype_model_reference, phenotype_model_effects) %>%
  mutate(
    phenotype = factor(phenotype, levels = phenotype_levels),
    model_outcome = factor(
      model_outcome,
      levels = c("Death or hospice", "Post-72h death or hospice", "ICU length of stay")
    ),
    estimand_label = paste0(as.character(model_outcome), "\n", estimand),
    estimand_label = factor(
      estimand_label,
      levels = c(
        "Death or hospice\nAdjusted odds ratio",
        "Post-72h death or hospice\nAdjusted hazard ratio",
        "ICU length of stay\nAdjusted mean ratio"
      )
    ),
    label_x = if_else(estimate_ci == "Reference", 1.05, pooled_conf_high_random * 1.06)
  )

readr::write_csv(
  phenotype_model_plot_df,
  file.path(latest_pooled_dir, "publication_table_phenotype_adjusted_model_effects.csv")
)

p_phenotype_model_forest <- ggplot(
  phenotype_model_plot_df,
  aes(
    x = pooled_effect_random,
    y = phenotype,
    xmin = pooled_conf_low_random,
    xmax = pooled_conf_high_random,
    color = phenotype
  )
) +
  geom_vline(xintercept = 1, linetype = 2, color = "gray45") +
  geom_errorbar(width = 0.14, linewidth = 0.9, orientation = "y") +
  geom_point(size = 4.2) +
  geom_text(
    aes(x = label_x, label = estimate_ci),
    hjust = 0,
    size = 4.1,
    color = "gray15",
    show.legend = FALSE
  ) +
  facet_wrap(~ estimand_label, ncol = 1, scales = "free_x") +
  scale_y_discrete(limits = rev(phenotype_levels)) +
  scale_color_manual(values = phenotype_palette, guide = "none") +
  scale_x_log10(
    breaks = c(0.8, 1, 1.5, 2, 3, 5, 8, 12, 16, 20),
    labels = c("0.8", "1.0", "1.5", "2", "3", "5", "8", "12", "16", "20"),
    expand = expansion(mult = c(0.03, 0.22))
  ) +
  labs(
    title = "Adjusted Outcome Associations by Respiratory Trajectory Phenotype",
    x = "Adjusted ratio vs Resolved (random-effects, log scale)",
    y = NULL
  ) +
  theme_pub(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 17, face = "bold"),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 15, face = "bold"),
    plot.margin = margin(10, 65, 10, 10)
  )

save_plot(p_phenotype_model_forest, "forest_phenotype_adjusted_outcome_models.png", w = 12.5, h = 11.5)

cat("Saved phenotype outcome figures to: ", latest_pooled_dir, "\n", sep = "")
