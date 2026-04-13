suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

repo_dir <- normalizePath(".", mustWork = TRUE)

pooled_dirs <- list.dirs(file.path(repo_dir, "output"), recursive = FALSE, full.names = TRUE)
pooled_dirs <- pooled_dirs[grepl("pooled_", basename(pooled_dirs))]
stopifnot(length(pooled_dirs) > 0)
latest_pooled_dir <- pooled_dirs[order(file.info(pooled_dirs)$mtime, decreasing = TRUE)][1]

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

save_plot <- function(p, filename, w = 10, h = 7, dpi = 320) {
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

normalize_binary <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("1", "TRUE", "Yes", "yes") ~ "Yes",
    x %in% c("0", "FALSE", "No", "no") ~ "No",
    TRUE ~ x
  )
}

clean_device <- function(x) {
  x <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    x %in% c("UNK", "UNKNOWN", "NA", "") ~ "Unknown / missing",
    x == "OTHER" ~ "Other",
    x == "HIGH FLOW NC" ~ "High-flow NC",
    x == "NIPPV" ~ "NIV",
    x == "IMV" ~ "IMV",
    x == "ROOM AIR" ~ "Room air",
    x == "NASAL CANNULA" ~ "Nasal cannula",
    x == "FACE MASK" ~ "Face mask",
    x == "TRACH COLLAR" ~ "Trach collar",
    x == "CPAP" ~ "CPAP",
    TRUE ~ stringr::str_to_sentence(x)
  )
}

phenotype_palette <- c(
  "Room air / No ARF" = "#2E8B57",
  "Low-flow O2 / No ARF" = "#D4A017",
  "NIV / No ARF" = "#2B6CB0",
  "IMV / Predom. no ARF" = "#C05621",
  "IMV / Hypoxemic/mixed ARF" = "#B83280"
)

overall_cat <- readr::read_csv(file.path(latest_pooled_dir, "pooled_table1_overall_categorical.csv"), show_col_types = FALSE)
overall_cont <- readr::read_csv(file.path(latest_pooled_dir, "pooled_table1_overall_continuous.csv"), show_col_types = FALSE)
by_cluster_cat <- readr::read_csv(file.path(latest_pooled_dir, "pooled_table1_by_cluster_categorical.csv"), show_col_types = FALSE)
by_cluster_cont <- readr::read_csv(file.path(latest_pooled_dir, "pooled_table1_by_cluster_continuous.csv"), show_col_types = FALSE)
table1_site <- readr::read_csv(file.path(latest_pooled_dir, "manuscript_table1.csv"), show_col_types = FALSE)

# Overall epidemiology summary table
overall_epi <- bind_rows(
  overall_cont %>%
    filter(variable %in% c("age_years", "charlson_score", "sofa_total", "icu_los_hours", "imv_hours_72h", "vaso_hours_72h", "pm25_5y", "no2_5y")) %>%
    transmute(
      domain = "Continuous",
      metric = recode(
        variable,
        age_years = "Age, mean (SD)",
        charlson_score = "Charlson score, mean (SD)",
        sofa_total = "SOFA total, mean (SD)",
        icu_los_hours = "ICU LOS (hours), mean (SD)",
        imv_hours_72h = "IMV hours in first 72h, mean (SD)",
        vaso_hours_72h = "Vasopressor hours in first 72h, mean (SD)",
        pm25_5y = "PM2.5 5-year exposure, mean (SD)",
        no2_5y = "NO2 5-year exposure, mean (SD)"
      ),
      value = sprintf("%.1f (%.1f)", mean, sd)
    ),
  overall_cat %>%
    filter(variable %in% c("death_in_hosp", "hospice_discharge", "death_or_hospice", "any_imv_72h", "vaso_any_72h",
                           "advanced_cancer_any_poa", "metastatic_cancer_poa", "malignant_pleural_effusion_poa", "cachexia_poa"),
           level %in% c("1", "TRUE")) %>%
    transmute(
      domain = "Categorical",
      metric = recode(
        variable,
        death_in_hosp = "In-hospital death, n (%)",
        hospice_discharge = "Hospice discharge, n (%)",
        death_or_hospice = "Death or hospice, n (%)",
        any_imv_72h = "Any IMV in first 72h, n (%)",
        vaso_any_72h = "Any vasopressor in first 72h, n (%)",
        advanced_cancer_any_poa = "Advanced cancer, n (%)",
        metastatic_cancer_poa = "Metastatic cancer, n (%)",
        malignant_pleural_effusion_poa = "Malignant pleural effusion, n (%)",
        cachexia_poa = "Cachexia, n (%)"
      ),
      value = sprintf("%s (%.1f%%)", comma(as.numeric(n)), 100 * as.numeric(prop))
    )
) %>%
  arrange(domain, metric)

readr::write_csv(overall_epi, file.path(latest_pooled_dir, "publication_table_overall_epi_summary.csv"))

# Site-level supplement table from manuscript table 1
site_epi_keep <- c(
  "N",
  "Age, mean (SD)",
  "Charlson score, mean (SD)",
  "SOFA total, mean (SD)",
  "Advanced cancer, n (%)",
  "In-hospital death, n (%)",
  "Hospice discharge, n (%)",
  "Death or hospice, n (%)",
  "Any vasopressor in first 72h, n (%)"
)

site_epi_table <- table1_site %>%
  filter(label %in% site_epi_keep)

readr::write_csv(site_epi_table, file.path(latest_pooled_dir, "supplement_table_site_epidemiology.csv"))

# Initial respiratory support overall
device_overall <- overall_cat %>%
  filter(variable == "device_t0") %>%
  mutate(device = clean_device(level)) %>%
  group_by(device) %>%
  summarize(
    n = sum(as.numeric(n), na.rm = TRUE),
    denom = max(as.numeric(denom), na.rm = TRUE),
    prop = n / denom,
    .groups = "drop"
  ) %>%
  arrange(desc(prop))

readr::write_csv(device_overall, file.path(latest_pooled_dir, "publication_table_initial_support_distribution.csv"))

p_device_overall <- device_overall %>%
  mutate(device = forcats::fct_reorder(device, prop)) %>%
  ggplot(aes(x = prop, y = device)) +
  geom_col(fill = "#2B6CB0") +
  geom_text(aes(label = sprintf("%s (%.1f%%)", comma(n), 100 * prop)), hjust = -0.05, size = 3.3) +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Initial respiratory support at ICU admission",
    x = "Percent of pooled cohort",
    y = NULL
  ) +
  theme_pub()

save_plot(p_device_overall, "figure_initial_support_distribution.png", w = 9, h = 6.5)

# Overall burden figure
burden_metrics <- overall_cat %>%
  filter(
    variable %in% c("advanced_cancer_any_poa", "metastatic_cancer_poa", "malignant_pleural_effusion_poa", "cachexia_poa",
                    "any_imv_72h", "vaso_any_72h", "death_in_hosp", "hospice_discharge", "death_or_hospice"),
    level %in% c("1", "TRUE")
  ) %>%
  transmute(
    metric = recode(
      variable,
      advanced_cancer_any_poa = "Advanced cancer",
      metastatic_cancer_poa = "Metastatic cancer",
      malignant_pleural_effusion_poa = "Pleural effusion",
      cachexia_poa = "Cachexia",
      any_imv_72h = "Any IMV by 72h",
      vaso_any_72h = "Any vasopressor by 72h",
      death_in_hosp = "In-hospital death",
      hospice_discharge = "Hospice discharge",
      death_or_hospice = "Death or hospice"
    ),
    n = as.numeric(n),
    prop = as.numeric(prop)
  ) %>%
  arrange(prop)

p_burden <- burden_metrics %>%
  mutate(metric = forcats::fct_reorder(metric, prop)) %>%
  ggplot(aes(x = prop, y = metric)) +
  geom_col(fill = "#7C3AED") +
  geom_text(aes(label = sprintf("%s (%.1f%%)", comma(n), 100 * prop)), hjust = -0.05, size = 3.3) +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Critical care burden and cancer severity in the pooled cohort",
    x = "Percent of pooled cohort",
    y = NULL
  ) +
  theme_pub()

save_plot(p_burden, "figure_overall_critical_care_burden.png", w = 9, h = 7)

# Cancer severity by phenotype
cancer_by_cluster <- by_cluster_cat %>%
  filter(
    variable %in% c("advanced_cancer_any_poa", "metastatic_cancer_poa", "malignant_pleural_effusion_poa", "cachexia_poa"),
    level %in% c("1", "TRUE")
  ) %>%
  transmute(
    phenotype = stringr::str_remove(pooled_cluster_label, "^Cluster \\d+: "),
    phenotype = recode(
      phenotype,
      "Low-flow oxygen / No ARF" = "Low-flow O2 / No ARF",
      "Noninvasive ventilation / No ARF" = "NIV / No ARF",
      "Invasive ventilation / Predominantly no ARF" = "IMV / Predom. no ARF",
      "Invasive ventilation / Hypoxemic or mixed ARF" = "IMV / Hypoxemic/mixed ARF"
    ),
    marker = recode(
      variable,
      advanced_cancer_any_poa = "Advanced cancer",
      metastatic_cancer_poa = "Metastatic cancer",
      malignant_pleural_effusion_poa = "Pleural effusion",
      cachexia_poa = "Cachexia"
    ),
    prop = as.numeric(prop),
    n = as.numeric(n)
  )

readr::write_csv(cancer_by_cluster, file.path(latest_pooled_dir, "publication_table_cancer_burden_by_cluster.csv"))

p_cancer_cluster <- ggplot(cancer_by_cluster, aes(x = prop, y = phenotype, color = marker)) +
  geom_point(size = 3) +
  geom_line(aes(group = marker), linewidth = 0.9, alpha = 0.8) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Cancer severity markers across pooled trajectory phenotypes",
    x = "Percent within phenotype",
    y = NULL,
    color = NULL
  ) +
  theme_pub() +
  theme(legend.position = "bottom")

save_plot(p_cancer_cluster, "figure_cancer_burden_by_cluster.png", w = 11, h = 6.5)

# Start support to pooled phenotype alluvial
alluvial_src <- by_cluster_cat %>%
  filter(variable == "device_t0") %>%
  transmute(
    start_support = clean_device(level),
    phenotype = stringr::str_remove(pooled_cluster_label, "^Cluster \\d+: "),
    phenotype = recode(
      phenotype,
      "Low-flow oxygen / No ARF" = "Low-flow O2 / No ARF",
      "Noninvasive ventilation / No ARF" = "NIV / No ARF",
      "Invasive ventilation / Predominantly no ARF" = "IMV / Predom. no ARF",
      "Invasive ventilation / Hypoxemic or mixed ARF" = "IMV / Hypoxemic/mixed ARF"
    ),
    n = as.numeric(n)
  ) %>%
  group_by(start_support, phenotype) %>%
  summarize(n = sum(n, na.rm = TRUE), .groups = "drop") %>%
  filter(n > 0)

readr::write_csv(alluvial_src, file.path(latest_pooled_dir, "publication_table_start_support_to_phenotype.csv"))

left_levels <- c("Room air", "Nasal cannula", "Face mask", "High-flow NC", "CPAP", "NIV", "IMV", "Trach collar", "Other", "Unknown / missing")
right_levels <- names(phenotype_palette)

flow_df <- alluvial_src %>%
  mutate(
    start_support = factor(start_support, levels = left_levels),
    phenotype = factor(phenotype, levels = right_levels)
  ) %>%
  arrange(start_support, phenotype)

left_pos <- flow_df %>%
  group_by(start_support) %>%
  summarize(total = sum(n), .groups = "drop") %>%
  arrange(start_support) %>%
  mutate(
    ymin = lag(cumsum(total), default = 0),
    ymax = cumsum(total)
  )

right_pos <- flow_df %>%
  group_by(phenotype) %>%
  summarize(total = sum(n), .groups = "drop") %>%
  arrange(phenotype) %>%
  mutate(
    ymin = lag(cumsum(total), default = 0),
    ymax = cumsum(total)
  )

flow_df <- flow_df %>%
  left_join(left_pos, by = "start_support", suffix = c("", "_left")) %>%
  rename(left_ymin = ymin, left_ymax = ymax) %>%
  left_join(right_pos, by = "phenotype", suffix = c("", "_right")) %>%
  rename(right_ymin = ymin, right_ymax = ymax) %>%
  group_by(start_support) %>%
  mutate(
    flow_left_ymin = left_ymin + cumsum(lag(n, default = 0)),
    flow_left_ymax = flow_left_ymin + n
  ) %>%
  ungroup() %>%
  arrange(phenotype, start_support) %>%
  group_by(phenotype) %>%
  mutate(
    flow_right_ymin = right_ymin + cumsum(lag(n, default = 0)),
    flow_right_ymax = flow_right_ymin + n
  ) %>%
  ungroup() %>%
  mutate(flow_id = row_number())

t_vals <- seq(0, 1, length.out = 60)
smoothstep <- 3 * t_vals^2 - 2 * t_vals^3

polys <- purrr::map_dfr(seq_len(nrow(flow_df)), function(i) {
  r <- flow_df[i, ]
  upper <- r$flow_left_ymax + (r$flow_right_ymax - r$flow_left_ymax) * smoothstep
  lower <- r$flow_left_ymin + (r$flow_right_ymin - r$flow_left_ymin) * smoothstep
  tibble(
    flow_id = r$flow_id,
    phenotype = as.character(r$phenotype),
    x = c(seq(1.02, 1.98, length.out = length(t_vals)),
          rev(seq(1.02, 1.98, length.out = length(t_vals)))),
    y = c(upper, rev(lower))
  )
})

left_rect <- left_pos %>%
  transmute(
    xmin = 0.92, xmax = 1.0,
    ymin, ymax,
    label = as.character(start_support)
  )

right_rect <- right_pos %>%
  transmute(
    xmin = 2.0, xmax = 2.08,
    ymin, ymax,
    phenotype = as.character(phenotype),
    label = as.character(phenotype)
  )

p_alluvial <- ggplot() +
  geom_polygon(data = polys, aes(x = x, y = y, group = flow_id, fill = phenotype), alpha = 0.78, color = NA) +
  geom_rect(data = left_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "grey20", color = NA) +
  geom_rect(data = right_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = phenotype), color = NA) +
  geom_text(data = left_rect, aes(x = 0.9, y = (ymin + ymax) / 2, label = label), hjust = 1, size = 3.1) +
  geom_text(data = right_rect, aes(x = 2.1, y = (ymin + ymax) / 2, label = label), hjust = 0, size = 3.1, fontface = "bold") +
  annotate("text", x = 0.78, y = max(right_rect$ymax) * 1.02, label = "Initial support", fontface = "bold", size = 4) +
  annotate("text", x = 2.24, y = max(right_rect$ymax) * 1.02, label = "Pooled phenotype", fontface = "bold", size = 4) +
  scale_fill_manual(values = phenotype_palette, drop = FALSE) +
  scale_y_continuous(labels = comma) +
  coord_cartesian(xlim = c(0.45, 2.55), clip = "off") +
  labs(
    title = "Initial respiratory support and pooled trajectory phenotype",
    x = NULL,
    y = "Hospitalizations",
    fill = "Pooled phenotype"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(10, 130, 10, 130)
  )

save_plot(p_alluvial, "figure_start_support_to_phenotype_alluvial.png", w = 14, h = 9)

cat("Saved lung cancer ICU epidemiology outputs to: ", latest_pooled_dir, "\n", sep = "")
