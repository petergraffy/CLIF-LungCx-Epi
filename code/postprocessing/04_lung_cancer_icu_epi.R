suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

repo_dir <- normalizePath(".", mustWork = TRUE)
site_dir <- file.path(repo_dir, "sites")

pooled_dirs <- list.dirs(file.path(repo_dir, "output"), recursive = FALSE, full.names = TRUE)
pooled_dirs <- pooled_dirs[grepl("pooled_", basename(pooled_dirs))]
stopifnot(length(pooled_dirs) > 0)
latest_pooled_dir <- pooled_dirs[order(basename(pooled_dirs), decreasing = TRUE)][1]

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
    x %in% c("UNK", "UNKNOWN", "NA", "", "NAN") ~ "Other",
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

prototype_tbl <- readr::read_csv(file.path(latest_pooled_dir, "pooled_cluster_prototypes.csv"), show_col_types = FALSE) %>%
  transmute(
    pooled_cluster = as.integer(pooled_cluster),
    phenotype = phenotype_label
  ) %>%
  arrange(pooled_cluster)

phenotype_levels <- c("Minimal ARF", "Resolving ARF", "Stable ARF", "Persistent ARF")
phenotype_levels <- phenotype_levels[phenotype_levels %in% prototype_tbl$phenotype]

phenotype_palette_master <- c(
  "Minimal ARF" = "#2E8B57",
  "Stable ARF" = "#D4A017",
  "Resolving ARF" = "#2B6CB0",
  "Persistent ARF" = "#C05621"
)
phenotype_palette <- phenotype_palette_master[phenotype_levels]

start_support_palette <- c(
  "Room air" = "#5B8E7D",
  "Nasal cannula" = "#A6C36F",
  "Face mask" = "#F2C14E",
  "High-flow NC" = "#F78154",
  "CPAP" = "#7D5BA6",
  "NIV" = "#8E63CE",
  "IMV" = "#C1666B",
  "Trach collar" = "#7A9E9F",
  "Other" = "#6C757D"
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
  mutate(device = coalesce(clean_device(level), "Other")) %>%
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
    marker = recode(
      variable,
      advanced_cancer_any_poa = "Advanced cancer",
      metastatic_cancer_poa = "Metastatic cancer",
      malignant_pleural_effusion_poa = "Pleural effusion",
      cachexia_poa = "Cachexia"
    ),
    prop = as.numeric(prop),
    n = as.numeric(n)
  ) %>%
  mutate(phenotype = factor(phenotype, levels = phenotype_levels))

readr::write_csv(cancer_by_cluster, file.path(latest_pooled_dir, "publication_table_cancer_burden_by_cluster.csv"))

p_cancer_cluster <- ggplot(cancer_by_cluster, aes(x = marker, y = prop, fill = phenotype)) +
  geom_col(position = position_dodge(width = 0.76), width = 0.68, color = "white") +
  geom_text(
    aes(label = sprintf("%.1f%%", 100 * prop)),
    position = position_dodge(width = 0.76),
    vjust = -0.25,
    size = 2.9
  ) +
  scale_fill_manual(values = phenotype_palette, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Cancer severity markers across pooled trajectory phenotypes",
    x = NULL,
    y = "Percent within phenotype",
    fill = NULL
  ) +
  theme_pub() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    axis.ticks.y = element_line(color = "black")
  )

save_plot(p_cancer_cluster, "figure_cancer_burden_by_cluster.png", w = 11, h = 6.5)

# Start support to pooled phenotype alluvial
alluvial_src <- by_cluster_cat %>%
  filter(variable == "device_t0") %>%
  transmute(
    start_support = coalesce(clean_device(level), "Other"),
    phenotype = stringr::str_remove(pooled_cluster_label, "^Cluster \\d+: "),
    n = as.numeric(n)
  ) %>%
  group_by(start_support, phenotype) %>%
  summarize(n = sum(n, na.rm = TRUE), .groups = "drop") %>%
  filter(n > 0)

readr::write_csv(alluvial_src, file.path(latest_pooled_dir, "publication_table_start_support_to_phenotype.csv"))

left_levels <- c(
  "Room air",
  "Nasal cannula",
  "Face mask",
  "High-flow NC",
  "CPAP",
  "NIV",
  "IMV",
  "Trach collar",
  "Other"
)
right_levels <- phenotype_levels

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
    start_support = as.character(start_support),
    n_total = ymax - ymin,
    label = paste0(start_support, ", n=", comma(n_total))
  )

right_rect <- right_pos %>%
  transmute(
    xmin = 2.0, xmax = 2.08,
    ymin, ymax,
    phenotype = as.character(phenotype),
    n_total = ymax - ymin,
    label = paste0(phenotype, ", n=", comma(n_total))
  )

p_alluvial <- ggplot() +
  geom_polygon(data = polys, aes(x = x, y = y, group = flow_id, fill = phenotype), alpha = 0.78, color = NA) +
  geom_rect(data = left_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = start_support), color = "white", linewidth = 0.45) +
  geom_rect(data = right_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = phenotype), color = NA) +
  geom_text(data = left_rect, aes(x = 0.9, y = (ymin + ymax) / 2, label = label), hjust = 1, size = 3.1) +
  geom_text(data = right_rect, aes(x = 2.1, y = (ymin + ymax) / 2, label = label), hjust = 0, size = 3.1, fontface = "bold") +
  annotate("text", x = 0.78, y = max(right_rect$ymax) * 1.02, label = "Initial support", fontface = "bold", size = 4) +
  annotate("text", x = 2.24, y = max(right_rect$ymax) * 1.02, label = "Pooled phenotype", fontface = "bold", size = 4) +
  scale_fill_manual(values = c(start_support_palette, phenotype_palette), drop = FALSE) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0.45, 2.55), clip = "off") +
  labs(
    title = "Initial respiratory support and pooled trajectory phenotype",
    x = NULL,
    y = "Hospitalizations"
  ) +
  theme_pub() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.y = element_line(color = "black", linewidth = 0.6),
    axis.ticks.y = element_line(color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 130, 10, 130)
  )

save_plot(p_alluvial, "figure_start_support_to_phenotype_alluvial.png", w = 14, h = 9)

# Combined site-level admission year trends
admit_year_paths <- list.files(site_dir, pattern = "^epi_admission_year_summary.*\\.csv$", recursive = TRUE, full.names = TRUE)

clean_site_name <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("emory", "Emory") ~ "Emory",
    x %in% c("Hopkins") ~ "Hopkins",
    x %in% c("U_Michigan", "Michigan") ~ "Michigan",
    x %in% c("nu", "NU") ~ "NU",
    x %in% c("RUSH", "Rush") ~ "Rush",
    x %in% c("UCMC") ~ "UCMC",
    x %in% c("UCSF") ~ "UCSF",
    x %in% c("umn", "UMN") ~ "UMN",
    TRUE ~ x
  )
}

if (length(admit_year_paths) > 0) {
  site_year <- purrr::map_dfr(admit_year_paths, ~ readr::read_csv(.x, show_col_types = FALSE)) %>%
    mutate(site_name = clean_site_name(site_name)) %>%
    arrange(site_name, admit_year)

  # Minnesota 2024 appears to be a truncated partial-year extract (n=2), so carry forward 2023 values for plotting.
  site_year <- site_year %>%
    group_by(site_name) %>%
    mutate(
      across(
        c(n, death_or_hospice_pct, any_imv_72h_pct, vaso_any_72h_pct, advanced_cancer_pct, mean_sofa_total),
        ~ ifelse(site_name == "UMN" & admit_year == 2024 & n < 10, dplyr::lag(.x), .x)
      )
    ) %>%
    ungroup()

  pooled_year <- site_year %>%
    group_by(admit_year) %>%
    group_modify(~ {
      w <- .x$n
      tibble(
        site_name = "Pooled",
        n = mean(w, na.rm = TRUE),
        death_or_hospice_pct = weighted.mean(.x$death_or_hospice_pct, w = w, na.rm = TRUE),
        any_imv_72h_pct = weighted.mean(.x$any_imv_72h_pct, w = w, na.rm = TRUE),
        vaso_any_72h_pct = weighted.mean(.x$vaso_any_72h_pct, w = w, na.rm = TRUE),
        advanced_cancer_pct = weighted.mean(.x$advanced_cancer_pct, w = w, na.rm = TRUE),
        mean_sofa_total = weighted.mean(.x$mean_sofa_total, w = w, na.rm = TRUE)
      )
    }) %>%
    ungroup()

  site_year_plot <- bind_rows(site_year, pooled_year) %>%
    pivot_longer(
      cols = c(n, death_or_hospice_pct, any_imv_72h_pct, vaso_any_72h_pct, mean_sofa_total),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric = factor(
        metric,
        levels = c("n", "death_or_hospice_pct", "any_imv_72h_pct", "vaso_any_72h_pct", "mean_sofa_total"),
        labels = c(
          "Hospitalizations",
          "Death or hospice discharge (%)",
          "Any IMV in first 72h (%)",
          "Any vasopressor in first 72h (%)",
          "Mean SOFA total"
        )
      ),
      value = ifelse(grepl("\\(\\%\\)$", as.character(metric)), 100 * value, value),
      site_name = factor(site_name, levels = c(sort(unique(setdiff(site_name, "Pooled"))), "Pooled"))
    )

  readr::write_csv(site_year_plot, file.path(latest_pooled_dir, "publication_table_site_admission_year_trends.csv"))

  p_admission_year <- ggplot(site_year_plot, aes(x = admit_year, y = value, group = site_name, color = site_name)) +
    geom_line(aes(linewidth = site_name == "Pooled"), alpha = 0.9) +
    geom_point(aes(size = site_name == "Pooled")) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    scale_linewidth_manual(values = c(`FALSE` = 0.8, `TRUE` = 1.6), guide = "none") +
    scale_size_manual(values = c(`FALSE` = 1.6, `TRUE` = 2.8), guide = "none") +
    scale_x_continuous(breaks = sort(unique(site_year_plot$admit_year))) +
    scale_color_manual(
      values = c(
        "Emory" = "#E76F51",
        "Hopkins" = "#E9C46A",
        "Michigan" = "#8D99AE",
        "NU" = "#2A9D8F",
        "Penn" = "#5C4D7D",
        "Rush" = "#F4A261",
        "UCMC" = "#3A86FF",
        "UCSF" = "#C77DFF",
        "UMN" = "#43AA8B",
        "Pooled" = "#111111"
      )
    ) +
    scale_y_continuous(labels = label_number(accuracy = 1)) +
    labs(
      title = "Admission year trends across sites",
      x = "Admission year",
      y = NULL,
      color = "Site"
    ) +
    theme_pub() +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.6),
      axis.line.y = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black")
    )

  save_plot(p_admission_year, "figure_admission_year_site_trends.png", w = 12, h = 9)
}

cat("Saved lung cancer ICU epidemiology outputs to: ", latest_pooled_dir, "\n", sep = "")
