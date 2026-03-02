# ================================================================================================
# CLIF | Lung Cancer ICU Cohort — Hourly Trajectory Phenotypes + Air Pollution Association
#
# What this script does:
#   1) Builds an hourly (0–72h post-ICU t0) respiratory support state panel from CLIF respiratory_support
#   2) Derives hourly vasoactive use from CLIF medication_admin_continuous (vasoactives; administered)
#   3) Derives hourly physiologic ARF subtypes (hypoxemic / hypercapnic / mixed / none) using:
#        - SpO2 paired to nearest FiO2 (room-air hypoxemia)
#        - PaO2 paired to nearest FiO2 (room-air hypoxemia + P/F <= 300)
#        - pCO2 paired to nearest pH (hypercapnic ARF)
#      with optional carry-forward smoothing for sparse ABGs
#   4) Constructs a composite hourly multi-organ failure (MOF) state:
#        RS_STATE | ARF_SUBTYPE | VASO
#   5) (Next steps) Performs TraMineR sequence clustering on composite MOF states
#   6) (Next steps) Models air pollution exposure → trajectory phenotype and → outcomes
#
# Assumptions:
#   - cohort_lung has hospitalization_id and t0 (POSIXct UTC)
#   - rs_raw (clif_respiratory_support) has recorded_dttm, device_category, fio2_set
#   - clif_tables has clif_vitals and clif_labs (or you load vitals/labs similarly)
#   - device_state_clif() and device_rank_clif() are defined to map device_category → ordered tiers
#   - analysis_ready already contains pm25/no2 windows and outcomes (death_or_hospice, icu_los_hours)
# ================================================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(TraMineR)
  library(cluster)
  library(data.table)
  library(nnet)
  library(ggplot2)
})

# ---------------------------
# Global params
# ---------------------------
H_MAX <- 72L
ROOM_AIR_FIO2 <- 0.21
JOIN_NEAR_H   <- 1   # hours for SpO2/PaO2 ↔ FiO2 nearest pairing
HYPERPAIR_H   <- 2   # hours for pCO2 ↔ pH nearest pairing
CARRY_H       <- 6L  # carry-forward window for hypox/hyper flags (ABG sparsity)

# ---------------------------
# Keys / anchor (one place)
# ---------------------------
keys <- cohort_lung %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    t0 = as.POSIXct(t0, tz="UTC")
  )

win72 <- keys %>%
  transmute(
    hospitalization_id,
    win_start = t0,
    win_end   = t0 + dhours(H_MAX)
  )

# ================================================================================================
# 1) Hourly Respiratory Support Panel (RS_STATE)
# ================================================================================================

rs_for_traj <- rs_raw %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm      = safe_posix(recorded_dttm),
    device_category    = as.character(device_category)
  ) %>%
  inner_join(keys, by="hospitalization_id") %>%
  mutate(
    dt_h = as.numeric(difftime(recorded_dttm, t0, units="hours")),
    h    = floor(dt_h),
    rs_state = device_state_clif(device_category),
    rank     = device_rank_clif(device_category)
  ) %>%
  filter(!is.na(h), h >= 0, h <= H_MAX, !is.na(rank))

hourly_rs <- rs_for_traj %>%
  arrange(hospitalization_id, h, desc(rank), recorded_dttm) %>%
  group_by(hospitalization_id, h) %>%
  slice(1) %>%
  ungroup() %>%
  select(hospitalization_id, h, rs_state)

rs_panel <- hourly_rs %>%
  group_by(hospitalization_id) %>%
  complete(h = 0:H_MAX) %>%
  arrange(hospitalization_id, h) %>%
  tidyr::fill(rs_state, .direction = "down") %>%
  tidyr::fill(rs_state, .direction = "up") %>%
  mutate(rs_state = coalesce(rs_state, "UNK")) %>%
  ungroup()

collapse_rs <- function(rs_state) {
  x <- toupper(trimws(as.character(rs_state)))
  
  case_when(
    x %in% c("ROOM AIR") ~ "ROOM AIR",
    
    x %in% c("NASAL CANNULA", "FACE MASK", "TRACH COLLAR") ~ "LOW_O2",
    
    x %in% c("HIGH FLOW NC", "CPAP", "NIPPV") ~ "NIV",
    
    x %in% c("IMV") ~ "IMV",
    
    # collapse UNKNOWN + OTHER into OTHER
    x %in% c("OTHER", "UNK", "", NA_character_) ~ "OTHER",
    
    TRUE ~ "OTHER"
  )
}

rs_panel <- rs_panel %>%
  mutate(rs_state_collapsed = collapse_rs(rs_state))

rs_panel %>% count(rs_state, rs_state_collapsed, sort = TRUE) %>% print(n = 50)

# ================================================================================================
# 2) Hourly Vasoactive Use (VASO)
#    Source: medication_admin_continuous where med_group == "vasoactives" and mar_action_group == "administered"
# ================================================================================================

med_cont <- get_min("medication_admin_continuous", c(
  "hospitalization_id","admin_dttm","med_group","mar_action_group"
)) %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    admin_dttm         = safe_posix(admin_dttm),
    med_group          = tolower(as.character(med_group)),
    mar_action_group   = tolower(as.character(mar_action_group))
  )

vaso_hourly <- med_cont %>%
  filter(
    med_group == "vasoactives",
    mar_action_group == "administered",
    !is.na(admin_dttm)
  ) %>%
  inner_join(keys, by="hospitalization_id") %>%
  mutate(
    dt_h = as.numeric(difftime(admin_dttm, t0, units="hours")),
    h    = floor(dt_h)
  ) %>%
  filter(!is.na(h), h >= 0, h <= H_MAX) %>%
  distinct(hospitalization_id, h) %>%
  mutate(vaso_h = 1L)

# ================================================================================================
# 3) Hourly ARF Subtype (ARF_HYPOX / ARF_HYPER / ARF_MIXED / NO_ARF)
#    Uses your physiologic definition, but time-resolved (0–72h window)
# ================================================================================================

# Inputs (long tables)
vitals <- clif_tables[["clif_vitals"]]
labs   <- clif_tables[["clif_labs"]]

# --- SpO2 observations (long vitals) ---
spo2_win <- vitals %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm      = safe_posix(recorded_dttm),
    vital_category     = tolower(as.character(vital_category)),
    vital_value        = as_num(vital_value)
  ) %>%
  filter(vital_category == "spo2") %>%
  inner_join(win72, by="hospitalization_id") %>%
  filter(recorded_dttm >= win_start, recorded_dttm <= win_end) %>%
  transmute(hospitalization_id, spo2_time = recorded_dttm, spo2 = vital_value)

# --- ABG labs (long labs) ---
labs_win <- labs %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    lab_result_dttm    = safe_posix(lab_result_dttm),
    lab_category       = tolower(as.character(lab_category)),
    lab_value_numeric  = as_num(lab_value_numeric)
  ) %>%
  filter(lab_category %in% c("po2_arterial","pco2_arterial","ph_arterial")) %>%
  inner_join(win72, by="hospitalization_id") %>%
  filter(lab_result_dttm >= win_start, lab_result_dttm <= win_end) %>%
  transmute(hospitalization_id, lab_time = lab_result_dttm, lab_category, val = lab_value_numeric)

# --- FiO2 set (from respiratory support) ---
fio2_win <- rs_raw %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm      = safe_posix(recorded_dttm),
    fio2_set           = as_num(fio2_set)
  ) %>%
  inner_join(win72, by="hospitalization_id") %>%
  filter(recorded_dttm >= win_start, recorded_dttm <= win_end) %>%
  transmute(hospitalization_id, fio2_time = recorded_dttm, fio2_set)

# --- Pairing via data.table nearest joins ---
setDT(spo2_win); setDT(fio2_win); setDT(labs_win)
setkey(spo2_win, hospitalization_id, spo2_time)
setkey(fio2_win, hospitalization_id, fio2_time)

# SpO2 ↔ FiO2 nearest within JOIN_NEAR_H
spo2_win[, spo2_time_keep := spo2_time]
spo2_fio2 <- fio2_win[
  spo2_win, roll="nearest", on=.(hospitalization_id, fio2_time = spo2_time), nomatch=0L
][
  , timediff_h := abs(as.numeric(difftime(spo2_time_keep, fio2_time, units="hours")))
][
  timediff_h <= as.numeric(JOIN_NEAR_H),
  .(hospitalization_id, t = spo2_time_keep, spo2, fio2_set,
    on_room_air = !is.na(fio2_set) & fio2_set <= ROOM_AIR_FIO2 + 1e-6)
]

# PaO2 ↔ FiO2 nearest within JOIN_NEAR_H
po2dt <- labs_win[lab_category == "po2_arterial",
                  .(hospitalization_id, po2_time = lab_time, po2 = val)]
setkey(po2dt, hospitalization_id, po2_time)
po2dt[, po2_time_keep := po2_time]
setkey(fio2_win, hospitalization_id, fio2_time)

po2_fio2 <- fio2_win[
  po2dt, roll="nearest", on=.(hospitalization_id, fio2_time = po2_time), nomatch=0L
][
  , timediff_h := abs(as.numeric(difftime(po2_time_keep, fio2_time, units="hours")))
][
  timediff_h <= as.numeric(JOIN_NEAR_H),
  .(hospitalization_id, t = po2_time_keep, po2, fio2_set,
    pf_ratio = fifelse(!is.na(fio2_set) & fio2_set > 0, po2 / fio2_set, as.numeric(NA)))
]

# pCO2 ↔ pH nearest within HYPERPAIR_H
pco2dt <- labs_win[lab_category == "pco2_arterial",
                   .(hospitalization_id, pco2_time = lab_time, pco2 = val)]
phdt   <- labs_win[lab_category == "ph_arterial",
                   .(hospitalization_id, ph_time = lab_time, ph = val)]
setkey(pco2dt, hospitalization_id, pco2_time)
setkey(phdt,   hospitalization_id, ph_time)

pco2dt[, pco2_time_keep := pco2_time]
phdt[,   ph_time_keep   := ph_time]

hyper_pairs <- phdt[
  pco2dt, roll="nearest", on=.(hospitalization_id, ph_time = pco2_time), nomatch=0L
][
  , timediff_h := abs(as.numeric(difftime(pco2_time_keep, ph_time, units="hours")))
][
  timediff_h <= as.numeric(HYPERPAIR_H),
  .(hospitalization_id, t = pco2_time_keep, pco2, ph,
    hyper_hit = (pco2 >= 45 & ph < 7.35))
]

# Helper: convert event timestamps to hour bins relative to t0
keys_dt <- as.data.table(keys)
setkey(keys_dt, hospitalization_id)

to_hour_bin <- function(dt, time_col) {
  dt[keys_dt, on="hospitalization_id"][,
                                       h := floor(as.numeric(difftime(get(time_col), t0, units="hours")))
  ][h >= 0 & h <= H_MAX]
}

spo2_fio2_h <- to_hour_bin(spo2_fio2, "t")[,
                                           .(hypox_roomair_spo2 = any(spo2 < 90 & on_room_air, na.rm=TRUE)),
                                           by=.(hospitalization_id, h)
]

po2_fio2_h <- to_hour_bin(po2_fio2, "t")[,
                                         .(
                                           hypox_roomair_po2 = any(po2 <= 60 & !is.na(fio2_set) & fio2_set <= ROOM_AIR_FIO2 + 1e-6, na.rm=TRUE),
                                           hypox_pf          = any(!is.na(pf_ratio) & pf_ratio <= 300, na.rm=TRUE)
                                         ),
                                         by=.(hospitalization_id, h)
]

hyper_h <- to_hour_bin(hyper_pairs, "t")[,
                                         .(hypercapnic = any(hyper_hit, na.rm=TRUE)),
                                         by=.(hospitalization_id, h)
]

arf_hourly <- Reduce(
  function(x, y) merge(x, y, by = c("hospitalization_id","h"), all = TRUE),
  list(spo2_fio2_h, po2_fio2_h, hyper_h)
) %>%
  as_tibble() %>%
  mutate(
    hypox_roomair_spo2 = coalesce(hypox_roomair_spo2, FALSE),
    hypox_roomair_po2  = coalesce(hypox_roomair_po2,  FALSE),
    hypox_pf           = coalesce(hypox_pf,           FALSE),
    hypercapnic        = coalesce(hypercapnic,        FALSE),
    hypox_h = hypox_roomair_spo2 | hypox_roomair_po2 | hypox_pf,
    hyper_h = hypercapnic
  ) %>%
  select(hospitalization_id, h, hypox_h, hyper_h)

carry_forward <- function(x, carry = 6L) {
  out <- rep(FALSE, length(x))
  idx <- which(x)
  for (i in idx) out[i:min(length(x), i + carry)] <- TRUE
  out
}

arf_hourly <- arf_hourly %>%
  group_by(hospitalization_id) %>%
  complete(h = 0:H_MAX, fill = list(hypox_h = FALSE, hyper_h = FALSE)) %>%
  arrange(h) %>%
  mutate(
    hypox_h = carry_forward(hypox_h, carry = CARRY_H),
    hyper_h = carry_forward(hyper_h, carry = CARRY_H),
    arf_subtype_h = case_when(
      hypox_h & hyper_h ~ "ARF_MIXED",
      hypox_h & !hyper_h ~ "ARF_HYPOX",
      !hypox_h & hyper_h ~ "ARF_HYPER",
      TRUE ~ "NO_ARF"
    )
  ) %>%
  ungroup() %>%
  select(hospitalization_id, h, arf_subtype_h)

# ================================================================================================
# 4) Composite MOF Panel: RS | ARF_SUBTYPE | VASO
# ================================================================================================

panel_mof <- rs_panel %>%
  left_join(arf_hourly, by=c("hospitalization_id","h")) %>%
  left_join(vaso_hourly, by=c("hospitalization_id","h")) %>%
  mutate(
    arf_subtype_h = coalesce(arf_subtype_h, "NO_ARF"),
    vaso_h = coalesce(vaso_h, 0L),
    state_mof = paste(
      rs_state,
      arf_subtype_h,
      ifelse(vaso_h == 1L, "VASO", "NoVASO"),
      sep="|"
    )
  )

# ================================================================================================
# 5) Sequence clustering on MOF composite states
# ================================================================================================
seq_wide_mof <- panel_mof %>%
  mutate(h_col = paste0("H", h)) %>%
  select(hospitalization_id, h_col, state_mof) %>%
  pivot_wider(names_from = h_col, values_from = state_mof)

seq_cols_mof <- seq_wide_mof %>% select(matches("^H\\d+$"))
alphabet_mof <- sort(unique(unlist(seq_cols_mof)))  # data-driven alphabet is safest here

seq_obj_mof <- TraMineR::seqdef(
  seq_cols_mof,
  alphabet = alphabet_mof,
  states   = alphabet_mof,
  labels   = alphabet_mof,
  xtstep   = 6
)

dist_mof <- TraMineR::seqdist(seq_obj_mof, method="OM", indel=1, sm="TRATE")
hc_mof <- hclust(as.dist(dist_mof), method="ward.D2")

K <- 8  # tune later
cluster_mof <- cutree(hc_mof, k=K)

traj_assign_mof <- tibble(
  hospitalization_id = seq_wide_mof$hospitalization_id,
  traj_cluster_mof = as.integer(cluster_mof)
)

save_csv(traj_assign_mof, "traj_cluster_assignments_mof")

png(file.path(out_dir, make_name("traj_mof_seqrplot_by_cluster", "png")), width=1600, height=1800, res=160)
TraMineR::seqrplot(seq_obj_mof, group = cluster_mof, diss = dist_mof,
                   main = "Representative MOF sequences by cluster (RS | ARF subtype | VASO)")
dev.off()


pick_k_silhouette <- function(dist_mat, k_min = 2, k_max = 12, hclust_method = "ward.D2") {
  stopifnot(is.matrix(dist_mat) || inherits(dist_mat, "dist"))
  d <- if (inherits(dist_mat, "dist")) dist_mat else as.dist(dist_mat)
  
  hc <- hclust(d, method = hclust_method)
  
  res <- lapply(k_min:k_max, function(k) {
    cl <- cutree(hc, k = k)
    sil <- silhouette(cl, d)
    tibble(
      k = k,
      avg_sil_width = mean(sil[, "sil_width"], na.rm = TRUE),
      min_sil_width = min(sil[, "sil_width"], na.rm = TRUE)
    )
  }) %>% bind_rows()
  
  p <- ggplot(res, aes(x = k, y = avg_sil_width)) +
    geom_line() + geom_point() +
    labs(title = "Average silhouette width vs K", x = "K (clusters)", y = "Avg silhouette width")
  
  list(table = res, plot = p, hc = hc)
}

# --- Example usage ---
# dist_mof <- TraMineR::seqdist(seq_obj_mof, method="OM", indel=1, sm="TRATE")
sil_out <- pick_k_silhouette(dist_mof, k_min = 2, k_max = 12)

print(sil_out$table)
print(sil_out$plot)

# "Best" K by silhouette
best_k_sil <- sil_out$table %>% slice_max(avg_sil_width, n = 1) %>% pull(k)
best_k_sil

panel_ra <- rs_panel %>%
  left_join(
    arf_hourly %>% select(hospitalization_id, h, arf_subtype_h),
    by = c("hospitalization_id","h")
  ) %>%
  mutate(
    arf_subtype_h = coalesce(arf_subtype_h, "NO_ARF"),
    rs_state_collapsed = coalesce(rs_state_collapsed, "OTHER"),
    state_ra = paste(rs_state_collapsed, arf_subtype_h, sep="|")
  )

seq_wide_ra <- panel_ra %>%
  mutate(h_col = paste0("H", h)) %>%
  select(hospitalization_id, h_col, state_ra) %>%
  pivot_wider(names_from = h_col, values_from = state_ra)

seq_cols_ra <- seq_wide_ra %>% select(matches("^H\\d+$"))
alphabet_ra <- sort(unique(unlist(seq_cols_ra)))

seq_obj_ra <- TraMineR::seqdef(
  seq_cols_ra,
  alphabet = alphabet_ra,
  states   = alphabet_ra,
  labels   = alphabet_ra,
  xtstep   = 6
)

dist_ra <- TraMineR::seqdist(seq_obj_ra, method="OM", indel=1, sm="TRATE")
hc_ra <- hclust(as.dist(dist_ra), method="ward.D2")

sil_out <- pick_k_silhouette(dist_ra, k_min = 2, k_max = 10)
print(sil_out$table)
print(sil_out$plot)

K <- 5   # or chosen from silhouette

cluster_ra <- cutree(hc_ra, k = K)

traj_assign_ra <- tibble(
  hospitalization_id = seq_wide_ra$hospitalization_id,
  traj_cluster_ra = as.integer(cluster_ra)
)

analysis_ready <- analysis_ready %>%
  left_join(traj_assign_ra, by="hospitalization_id") %>%
  mutate(traj_cluster_ra = factor(traj_cluster_ra))

png("traj_ra_seqrplot.png", width=2600, height=2800, res=160)
TraMineR::seqrplot(
  seq_obj_ra,
  group = cluster_ra,
  diss  = dist_ra,
  main = "Representative sequences (RS | ARF subtype)"
)
dev.off()


# Extract states
alphabet_ra <- sort(unique(unlist(seq_cols_ra)))

# Helper to assign color by RS tier + ARF subtype
state_colors <- function(state) {
  
  parts <- strsplit(state, "\\|")[[1]]
  rs  <- parts[1]
  arf <- parts[2]
  
  # Base color by RS tier
  base <- case_when(
    rs == "ROOM AIR" ~ "#2ca25f",   # green
    rs == "LOW_O2"   ~ "#e6ab02",   # gold
    rs == "NIV"      ~ "#3182bd",   # blue
    rs == "IMV"      ~ "#de2d26",   # red
    rs == "OTHER"    ~ "#636363",   # gray
    TRUE ~ "#969696"
  )
  
  # Modify brightness by ARF subtype
  adjust <- function(col, factor) {
    rgb_col <- col2rgb(col) / 255
    rgb_adj <- pmin(1, rgb_col * factor)
    rgb(rgb_adj[1], rgb_adj[2], rgb_adj[3])
  }
  
  case_when(
    arf == "NO_ARF"    ~ adjust(base, 1.2),
    arf == "ARF_HYPOX" ~ adjust(base, 1.0),
    arf == "ARF_HYPER" ~ adjust(base, 0.8),
    arf == "ARF_MIXED" ~ adjust(base, 0.6),
    TRUE ~ base
  )
}

cpal_ra <- sapply(alphabet_ra, state_colors)

seq_obj_ra <- TraMineR::seqdef(
  seq_cols_ra,
  alphabet = alphabet_ra,
  states   = alphabet_ra,
  labels   = alphabet_ra,
  cpal     = cpal_ra,
  xtstep   = 6
)

png("traj_ra_seqrplot_discrete.png", width=1600, height=2800, res=160)
TraMineR::seqrplot(
  seq_obj_ra,
  group = cluster_ra,
  diss  = dist_ra,
  main = "Representative sequences (Collapsed RS | ARF subtype)"
)
dev.off()










# ================================================================================================
# 5) Association to air pollution
# ================================================================================================


analysis_ready <- analysis_ready %>%
  mutate(hospitalization_id = as.character(hospitalization_id)) %>%
  left_join(traj_assign_mof, by="hospitalization_id")

analysis_ready %>% summarize(n=n(), missing=sum(is.na(traj_cluster_mof)))

analysis_ready <- analysis_ready %>%
  mutate(
    traj_cluster_mof = factor(traj_cluster_mof),
    pm25_5y_z = as.numeric(scale(pm25_5y)),
    no2_5y_z  = as.numeric(scale(no2_5y))
  )

m_traj_mof <- nnet::multinom(
  traj_cluster_mof ~ pm25_5y_z + no2_5y_z + age_years + sex_category + race_category,
  data = analysis_ready
)
summary(m_traj_mof)

m_mort_mof <- glm(
  death_or_hospice ~ pm25_5y_z + no2_5y_z + traj_cluster_mof +
    age_years + sex_category + race_category,
  data = analysis_ready,
  family = binomial()
)
summary(m_mort_mof)
exp(cbind(OR = coef(m_mort_mof), confint(m_mort_mof)))

m_int_pm <- glm(
  death_or_hospice ~ pm25_5y_z * traj_cluster_mof + no2_5y_z +
    age_years + sex_category + race_category,
  data = analysis_ready,
  family = binomial()
)

m_int_no2 <- glm(
  death_or_hospice ~ no2_5y_z * traj_cluster_mof + pm25_5y_z +
    age_years + sex_category + race_category,
  data = analysis_ready,
  family = binomial()
)
summary(m_int_pm)
summary(m_int_no2)

m_los_mof <- glm(
  icu_los_hours ~ pm25_5y_z + no2_5y_z + traj_cluster_mof +
    age_years + sex_category + race_category,
  data = analysis_ready,
  family = quasipoisson(link="log")
)
summary(m_los_mof)










analysis_ready %>% summarize(n=n(), missing=sum(is.na(traj_cluster_ra)))

analysis_ready <- analysis_ready %>%
  mutate(
    traj_cluster_ra = factor(traj_cluster_ra),
    pm25_5y_z = as.numeric(scale(pm25_5y)),
    no2_5y_z  = as.numeric(scale(no2_5y))
  )

m_traj_mof <- nnet::multinom(
  traj_cluster_ra ~ pm25_5y_z + no2_5y_z + age_years + sex_category + race_category,
  data = analysis_ready
)
summary(m_traj_mof)

m_mort_mof <- glm(
  death_or_hospice ~ pm25_5y_z + no2_5y_z + traj_cluster_ra +
    age_years + sex_category + race_category,
  data = analysis_ready,
  family = binomial()
)
summary(m_mort_mof)
exp(cbind(OR = coef(m_mort_mof), confint(m_mort_mof)))

m_int_pm <- glm(
  death_or_hospice ~ pm25_5y_z * traj_cluster_ra + no2_5y_z +
    age_years + sex_category + race_category,
  data = analysis_ready,
  family = binomial()
)

m_int_no2 <- glm(
  death_or_hospice ~ no2_5y_z * traj_cluster_ra + pm25_5y_z +
    age_years + sex_category + race_category,
  data = analysis_ready,
  family = binomial()
)
summary(m_int_pm)
summary(m_int_no2)

m_los_mof <- glm(
  icu_los_hours ~ pm25_5y_z + no2_5y_z + traj_cluster_ra +
    age_years + sex_category + race_category,
  data = analysis_ready,
  family = quasipoisson(link="log")
)
summary(m_los_mof)

car::vif(m_mort_mof)

analysis_ready %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    mortality = mean(death_or_hospice),
    pm25 = mean(pm25_5y, na.rm = TRUE),
    no2 = mean(no2_5y, na.rm = TRUE)
  )



# Build prediction grid
newdat_pm <- expand_grid(
  traj_cluster_ra = levels(analysis_ready$traj_cluster_ra),
  pm25_5y_z = seq(-2, 2, by = 0.1),
  no2_5y_z = 0,  # hold constant at mean
  age_years = mean(analysis_ready$age_years, na.rm = TRUE),
  sex_category = "Female",
  race_category = "White"
)

newdat_pm$pred_prob <- predict(
  m_int_pm,
  newdata = newdat_pm,
  type = "response"
)

p_pm <- ggplot(newdat_pm,
               aes(x = pm25_5y_z,
                   y = pred_prob,
                   color = traj_cluster_ra)) +
  geom_line(size = 1.2) +
  labs(
    x = "PM2.5 (SD units)",
    y = "Predicted Probability of Death",
    color = "Trajectory Cluster",
    title = "Cluster-Specific Mortality Risk by PM2.5 Exposure"
  ) +
  theme_minimal(base_size = 14)

p_pm

newdat_no2 <- expand_grid(
  traj_cluster_ra = levels(analysis_ready$traj_cluster_ra),
  no2_5y_z = seq(-2, 2, by = 0.1),
  pm25_5y_z = 0,
  age_years = mean(analysis_ready$age_years, na.rm = TRUE),
  sex_category = "Female",
  race_category = "White"
)

newdat_no2$pred_prob <- predict(
  m_int_no2,
  newdata = newdat_no2,
  type = "response"
)

p_no2 <- ggplot(newdat_no2,
                aes(x = no2_5y_z,
                    y = pred_prob,
                    color = traj_cluster_ra)) +
  geom_line(size = 1.2) +
  labs(
    x = "NO2 (SD units)",
    y = "Predicted Probability of Death",
    color = "Trajectory Cluster",
    title = "Cluster-Specific Mortality Risk by NO2 Exposure"
  ) +
  theme_minimal(base_size = 14)

p_no2


contrast_pm <- newdat_pm %>%
  filter(pm25_5y_z %in% c(-1, 1)) %>%
  mutate(pm25_label = ifelse(pm25_5y_z == -1, "low", "high")) %>%
  select(traj_cluster_ra, pm25_label, pred_prob) %>%
  pivot_wider(
    names_from = pm25_label,
    values_from = pred_prob
  ) %>%
  mutate(
    risk_diff = high - low,
    risk_ratio = high / low
  )

contrast_pm

contrast_no2 <- newdat_no2 %>%
  filter(no2_5y_z %in% c(-1, 1)) %>%
  mutate(no2_label = ifelse(no2_5y_z == -1, "low", "high")) %>%
  select(traj_cluster_ra, no2_label, pred_prob) %>%
  pivot_wider(
    names_from = no2_label,
    values_from = pred_prob
  ) %>%
  mutate(
    risk_diff = high - low,
    risk_ratio = high / low
  )

contrast_no2











