suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

repo_dir <- normalizePath(".", mustWork = TRUE)
pooled_dirs <- list.dirs(file.path(repo_dir, "output"), recursive = FALSE, full.names = TRUE)
pooled_dirs <- pooled_dirs[grepl("pooled_", basename(pooled_dirs))]
stopifnot(length(pooled_dirs) > 0)
latest_pooled_dir <- pooled_dirs[order(file.info(pooled_dirs)$mtime, decreasing = TRUE)][1]

mapping_path <- file.path(latest_pooled_dir, "pooled_cluster_mapping.csv")
stopifnot(file.exists(mapping_path))

out_dir <- latest_pooled_dir

theme_pub <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

save_plot <- function(p, filename, w = 13, h = 9, dpi = 320) {
  ggsave(
    filename = file.path(out_dir, filename),
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

mapping <- readr::read_csv(mapping_path, show_col_types = FALSE) %>%
  mutate(
    local_cluster = as.integer(local_cluster),
    local_severity_rank = as.integer(local_severity_rank),
    pooled_severity_rank = as.integer(pooled_severity_rank),
    n = as.numeric(n),
    phenotype = phenotype_label_short,
    site_name = factor(site_name, levels = unique(site_name)),
    local_label = paste0(site_name, " C", local_cluster)
  ) %>%
  arrange(site_name, local_severity_rank, local_cluster)

build_flow_polygons <- function(df, left_col, right_col, weight_col = "n", facet_col = NULL) {
  facet_quo <- rlang::enquo(facet_col)
  left_sym <- rlang::ensym(left_col)
  right_sym <- rlang::ensym(right_col)
  weight_sym <- rlang::ensym(weight_col)
  has_facet <- !rlang::quo_is_null(facet_quo)

  dat <- df %>%
    mutate(.facet = if (has_facet) as.character(!!facet_quo) else "All") %>%
    group_by(.facet, !!left_sym) %>%
    mutate(.left_total = sum(!!weight_sym, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(.facet, !!right_sym) %>%
    mutate(.right_total = sum(!!weight_sym, na.rm = TRUE)) %>%
    ungroup()

  left_levels <- dat %>%
    distinct(.facet, !!left_sym, .left_total) %>%
    arrange(.facet, !!left_sym) %>%
    group_by(.facet) %>%
    mutate(
      left_ymax = cumsum(.left_total),
      left_ymin = left_ymax - .left_total,
      left_y = (left_ymin + left_ymax) / 2
    ) %>%
    ungroup()

  right_levels <- dat %>%
    distinct(.facet, !!right_sym, .right_total) %>%
    arrange(.facet, !!right_sym) %>%
    group_by(.facet) %>%
    mutate(
      right_ymax = cumsum(.right_total),
      right_ymin = right_ymax - .right_total,
      right_y = (right_ymin + right_ymax) / 2
    ) %>%
    ungroup()

  flow_pos <- dat %>%
    left_join(left_levels, by = c(".facet", rlang::as_name(left_sym))) %>%
    left_join(right_levels, by = c(".facet", rlang::as_name(right_sym))) %>%
    arrange(.facet, !!left_sym, !!right_sym) %>%
    group_by(.facet, !!left_sym) %>%
    mutate(
      flow_left_ymin = left_ymin + cumsum(lag(!!weight_sym, default = 0)),
      flow_left_ymax = flow_left_ymin + !!weight_sym
    ) %>%
    ungroup() %>%
    arrange(.facet, !!right_sym, !!left_sym) %>%
    group_by(.facet, !!right_sym) %>%
    mutate(
      flow_right_ymin = right_ymin + cumsum(lag(!!weight_sym, default = 0)),
      flow_right_ymax = flow_right_ymin + !!weight_sym
    ) %>%
    ungroup() %>%
    mutate(flow_id = row_number())

  t_vals <- seq(0, 1, length.out = 60)
  smoothstep <- 3 * t_vals^2 - 2 * t_vals^3

  polys <- purrr::map_dfr(seq_len(nrow(flow_pos)), function(i) {
    r <- flow_pos[i, ]
    upper <- r$flow_left_ymax + (r$flow_right_ymax - r$flow_left_ymax) * smoothstep
    lower <- r$flow_left_ymin + (r$flow_right_ymin - r$flow_left_ymin) * smoothstep
    tibble(
      .facet = r$.facet,
      flow_id = r$flow_id,
      phenotype = r[[rlang::as_name(right_sym)]],
      x = c(seq(1.05, 1.95, length.out = length(t_vals)),
            rev(seq(1.05, 1.95, length.out = length(t_vals)))),
      y = c(upper, rev(lower))
    )
  })

  list(
    polygons = polys,
    left_strata = left_levels,
    right_strata = right_levels,
    flow_pos = flow_pos
  )
}

site_order <- mapping %>%
  distinct(site_name) %>%
  pull(site_name) %>%
  as.character()

left_order <- mapping %>%
  distinct(site_name, local_cluster, local_severity_rank, local_label) %>%
  arrange(factor(site_name, levels = site_order), local_severity_rank, local_cluster) %>%
  pull(local_label)

right_order <- mapping %>%
  distinct(pooled_severity_rank, phenotype) %>%
  arrange(pooled_severity_rank) %>%
  pull(phenotype)

alluvial_df <- mapping %>%
  mutate(
    local_label = factor(local_label, levels = left_order),
    phenotype = factor(phenotype, levels = right_order)
  )

flow_combined <- build_flow_polygons(alluvial_df, local_label, phenotype)

left_rect <- flow_combined$left_strata %>%
  transmute(
    .facet,
    xmin = 0.9,
    xmax = 1.0,
    ymin = left_ymin,
    ymax = left_ymax,
    label = as.character(local_label)
  )

right_rect <- flow_combined$right_strata %>%
  transmute(
    .facet,
    xmin = 2.0,
    xmax = 2.1,
    ymin = right_ymin,
    ymax = right_ymax,
    phenotype = as.character(phenotype),
    label = as.character(phenotype)
  )

p_combined <- ggplot() +
  geom_polygon(
    data = flow_combined$polygons,
    aes(x = x, y = y, group = flow_id, fill = phenotype),
    alpha = 0.75,
    color = NA
  ) +
  geom_rect(
    data = left_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "grey20",
    color = NA
  ) +
  geom_rect(
    data = right_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = phenotype),
    color = NA
  ) +
  geom_text(
    data = left_rect,
    aes(x = 0.88, y = (ymin + ymax) / 2, label = label),
    hjust = 1,
    size = 2.6,
    lineheight = 0.9
  ) +
  geom_text(
    data = right_rect,
    aes(x = 2.12, y = (ymin + ymax) / 2, label = label),
    hjust = 0,
    size = 3.2,
    lineheight = 0.95,
    fontface = "bold"
  ) +
  scale_fill_manual(values = phenotype_palette, drop = FALSE) +
  coord_cartesian(xlim = c(0.55, 2.55), clip = "off") +
  labs(
    title = "Mapping site-local trajectory clusters to pooled consensus phenotypes",
    x = NULL,
    y = "Hospitalizations",
    fill = "Pooled phenotype"
  ) +
  annotate("text", x = 0.75, y = max(right_rect$ymax) * 1.02, label = "Site-local clusters", fontface = "bold", size = 4) +
  annotate("text", x = 2.25, y = max(right_rect$ymax) * 1.02, label = "Consensus phenotypes", fontface = "bold", size = 4) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(labels = comma) +
  theme_pub() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(10, 130, 10, 170)
  )

save_plot(p_combined, "figure_alluvial_local_to_pooled.png", w = 16, h = 11)

facet_df <- mapping %>%
  mutate(
    local_label = factor(paste0("C", local_cluster), levels = unique(paste0("C", local_cluster))),
    phenotype = factor(phenotype, levels = right_order),
    site_name = factor(site_name, levels = site_order)
  )

flow_faceted <- build_flow_polygons(facet_df, local_label, phenotype, facet_col = site_name)

facet_left_rect <- flow_faceted$left_strata %>%
  transmute(
    .facet,
    xmin = 0.9,
    xmax = 1.0,
    ymin = left_ymin,
    ymax = left_ymax,
    label = as.character(local_label)
  )

facet_right_rect <- flow_faceted$right_strata %>%
  transmute(
    .facet,
    xmin = 2.0,
    xmax = 2.1,
    ymin = right_ymin,
    ymax = right_ymax,
    phenotype = as.character(phenotype),
    label = as.character(phenotype)
  )

p_faceted <- ggplot() +
  geom_polygon(
    data = flow_faceted$polygons,
    aes(x = x, y = y, group = flow_id, fill = phenotype),
    alpha = 0.78,
    color = NA
  ) +
  geom_rect(
    data = facet_left_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "grey20",
    color = NA
  ) +
  geom_rect(
    data = facet_right_rect,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = phenotype),
    color = NA
  ) +
  geom_text(
    data = facet_left_rect,
    aes(x = 0.88, y = (ymin + ymax) / 2, label = label),
    hjust = 1,
    size = 2.8
  ) +
  facet_wrap(~ .facet, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = phenotype_palette, drop = FALSE) +
  coord_cartesian(xlim = c(0.55, 2.35), clip = "off") +
  labs(
    title = "Site-specific local-to-pooled cluster mappings",
    x = NULL,
    y = "Hospitalizations",
    fill = "Pooled phenotype"
  ) +
  annotate("text", x = 0.75, y = Inf, label = "Local", fontface = "bold", vjust = 1.5, size = 3.4) +
  annotate("text", x = 2.18, y = Inf, label = "Pooled", fontface = "bold", vjust = 1.5, size = 3.4) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(labels = comma) +
  theme_pub() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 11),
    plot.margin = margin(10, 80, 10, 55)
  )

save_plot(p_faceted, "figure_alluvial_local_to_pooled_faceted.png", w = 15, h = 14)

readr::write_csv(
  alluvial_df %>%
    dplyr::select(site_name, local_cluster, local_severity_rank, pooled_cluster, pooled_severity_rank,
                  phenotype, n, death_pct, imv72_pct, sofa_total_mean, icu_los_median),
  file.path(out_dir, "alluvial_mapping_table.csv")
)

cat("Saved alluvial mockups to: ", out_dir, "\n", sep = "")
