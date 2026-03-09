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

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(ggplot2)
})

H_MAX <- 72L

# Ensure cluster is present
stopifnot("traj_cluster_ra" %in% names(analysis_ready))

# Build the per-hour panel with collapsed RS + ARF subtype and cluster assignment
sig_panel <- rs_panel %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    h = as.integer(h),
    rs = as.character(rs_state_collapsed)
  ) %>%
  left_join(
    arf_hourly %>%
      transmute(hospitalization_id = as.character(hospitalization_id),
                h = as.integer(h),
                arf = as.character(arf_subtype_h)),
    by = c("hospitalization_id","h")
  ) %>%
  left_join(
    analysis_ready %>% transmute(hospitalization_id = as.character(hospitalization_id),
                                 traj_cluster_ra),
    by = "hospitalization_id"
  ) %>%
  mutate(
    rs  = toupper(coalesce(rs, "OTHER")),
    arf = toupper(coalesce(arf, "NO_ARF")),
    traj_cluster_ra = forcats::fct_drop(traj_cluster_ra)
  ) %>%
  filter(!is.na(traj_cluster_ra), !is.na(h), h >= 0, h <= H_MAX)

# Quick QC
sig_panel %>% count(rs, sort = TRUE) %>% print(n = 50)
sig_panel %>% count(arf, sort = TRUE) %>% print(n = 50)

# RS signature: proportions by hour within each cluster
sig_rs <- sig_panel %>%
  count(traj_cluster_ra, h, rs) %>%
  group_by(traj_cluster_ra, h) %>%
  mutate(p = n / sum(n)) %>%
  ungroup()

# Set RS order (severity-ish)
rs_levels <- c("ROOM AIR", "LOW_O2", "NIV", "IMV", "OTHER")
sig_rs$rs <- factor(sig_rs$rs, levels = rs_levels)

p_rs <- ggplot(sig_rs, aes(x = h, y = p, fill = rs)) +
  geom_area() +
  facet_wrap(~ traj_cluster_ra, ncol = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Trajectory signatures by cluster: Respiratory support (collapsed tiers)",
    x = "Hour from ICU t0",
    y = "Within-cluster proportion"
  ) +
  theme_minimal(base_size = 14)

p_rs
save_plot(p_rs, "sig_rs_area_by_cluster", w = 12, h = 8)


# RS signature: proportions by hour within each cluster
sig_rs <- sig_panel %>%
  count(traj_cluster_ra, h, rs) %>%
  group_by(traj_cluster_ra, h) %>%
  mutate(p = n / sum(n)) %>%
  ungroup()

# Set RS order (severity-ish)
rs_levels <- c("ROOM AIR", "LOW_O2", "NIV", "IMV", "OTHER")
sig_rs$rs <- factor(sig_rs$rs, levels = rs_levels)

p_rs <- ggplot(sig_rs, aes(x = h, y = p, fill = rs)) +
  geom_area() +
  facet_wrap(~ traj_cluster_ra, ncol = 2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Trajectory signatures by cluster: Respiratory support (collapsed tiers)",
    x = "Hour from ICU t0",
    y = "Within-cluster proportion"
  ) +
  theme_minimal(base_size = 14)

p_rs
save_plot(p_rs, "sig_rs_area_by_cluster", w = 12, h = 8)

# patient-level first-hit times
first_hits <- sig_panel %>%
  group_by(hospitalization_id, traj_cluster_ra) %>%
  summarize(
    t_imv = suppressWarnings(min(h[rs == "IMV"], na.rm = TRUE)),
    t_arf = suppressWarnings(min(h[arf != "NO_ARF"], na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    t_imv = ifelse(is.infinite(t_imv), NA, t_imv),
    t_arf = ifelse(is.infinite(t_arf), NA, t_arf)
  )

# helper: cumulative incidence curve
cum_curve <- function(df, time_var, label) {
  df %>%
    filter(!is.na(.data[[time_var]])) %>%
    count(traj_cluster_ra, t = .data[[time_var]]) %>%
    group_by(traj_cluster_ra) %>%
    arrange(t) %>%
    mutate(cum = cumsum(n) / sum(n),
           outcome = label) %>%
    ungroup()
}

cc_imv <- cum_curve(first_hits, "t_imv", "IMV start")
cc_arf <- cum_curve(first_hits, "t_arf", "Any ARF onset")

cc <- bind_rows(cc_imv, cc_arf)

p_cc <- ggplot(cc, aes(x = t, y = cum, color = traj_cluster_ra)) +
  geom_step(linewidth = 1) +
  facet_wrap(~ outcome, ncol = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Timing signatures by cluster",
    x = "Hour from ICU t0",
    y = "Cumulative proportion with event",
    color = "Cluster"
  ) +
  theme_minimal(base_size = 14)

p_cc
save_plot(p_cc, "sig_timing_cuminc_by_cluster", w = 10, h = 8)


# Save previously created marginal-effect plots
save_plot(p_pm,  "marginal_pm25_by_cluster", w = 10, h = 7)
save_plot(p_no2, "marginal_no2_by_cluster",  w = 10, h = 7)

# Prediction grid: PM2.5
newdat_pm <- tidyr::expand_grid(
  traj_cluster_ra = levels(analysis_ready$traj_cluster_ra),
  pm25_5y_z = seq(-2, 2, by = 0.1),
  no2_5y_z = 0,
  age_years = mean(analysis_ready$age_years, na.rm = TRUE),
  sex_category = "Female",
  race_category = "White"
)
newdat_pm$pred_prob <- predict(m_int_pm, newdata = newdat_pm, type = "response")

p_pm <- ggplot(newdat_pm, aes(x = pm25_5y_z, y = pred_prob, color = traj_cluster_ra)) +
  geom_line(linewidth = 1.1) +
  labs(
    title = "Cluster-specific predicted mortality by PM2.5 (5y z-score)",
    x = "PM2.5 (SD units)",
    y = "Predicted probability of death/hospice",
    color = "Cluster"
  ) +
  theme_minimal(base_size = 14)

p_pm
save_plot(p_pm, "marginal_pm25_by_cluster", w = 10, h = 7)

# Prediction grid: NO2
newdat_no2 <- tidyr::expand_grid(
  traj_cluster_ra = levels(analysis_ready$traj_cluster_ra),
  no2_5y_z = seq(-2, 2, by = 0.1),
  pm25_5y_z = 0,
  age_years = mean(analysis_ready$age_years, na.rm = TRUE),
  sex_category = "Female",
  race_category = "White"
)
newdat_no2$pred_prob <- predict(m_int_no2, newdata = newdat_no2, type = "response")

p_no2 <- ggplot(newdat_no2, aes(x = no2_5y_z, y = pred_prob, color = traj_cluster_ra)) +
  geom_line(linewidth = 1.1) +
  labs(
    title = "Cluster-specific predicted mortality by NO2 (5y z-score)",
    x = "NO2 (SD units)",
    y = "Predicted probability of death/hospice",
    color = "Cluster"
  ) +
  theme_minimal(base_size = 14)

p_no2
save_plot(p_no2, "marginal_no2_by_cluster", w = 10, h = 7)

contrast_pm <- newdat_pm %>%
  filter(pm25_5y_z %in% c(-1, 1)) %>%
  mutate(pm25_label = ifelse(pm25_5y_z == -1, "low", "high")) %>%
  dplyr::select(traj_cluster_ra, pm25_label, pred_prob) %>%
  pivot_wider(names_from = pm25_label, values_from = pred_prob) %>%
  mutate(risk_diff = high - low, risk_ratio = high / low)

save_csv(contrast_pm, "contrast_pm25_low_vs_high_by_cluster")

contrast_no2 <- newdat_no2 %>%
  filter(no2_5y_z %in% c(-1, 1)) %>%
  mutate(no2_label = ifelse(no2_5y_z == -1, "low", "high")) %>%
  dplyr::select(traj_cluster_ra, no2_label, pred_prob) %>%
  pivot_wider(names_from = no2_label, values_from = pred_prob) %>%
  mutate(risk_diff = high - low, risk_ratio = high / low)

save_csv(contrast_no2, "contrast_no2_low_vs_high_by_cluster")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
})

# Prediction grid (hold covariates constant)
newdat_pm <- tidyr::expand_grid(
  traj_cluster_ra = levels(droplevels(analysis_ready$traj_cluster_ra)),
  pm25_5y_z = seq(-2, 2, by = 0.1),
  no2_5y_z = 0,
  age_years = mean(analysis_ready$age_years, na.rm = TRUE),
  sex_category = "Female",
  race_category = "White"
)

newdat_pm$pred <- predict(m_int_pm, newdata = newdat_pm, type = "response")

p_int_pm <- ggplot(newdat_pm, aes(x = pm25_5y_z, y = pred, color = traj_cluster_ra)) +
  geom_line(linewidth = 1.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Interaction: PM2.5 × Trajectory cluster (predicted mortality)",
    x = "PM2.5 (5y z-score)",
    y = "Predicted probability of death/hospice",
    color = "Cluster"
  ) +
  theme_minimal(base_size = 14)

p_int_pm
save_plot(p_int_pm, "interaction_pm25_by_cluster_predprob", w = 10, h = 7)

newdat_no2 <- tidyr::expand_grid(
  traj_cluster_ra = levels(droplevels(analysis_ready$traj_cluster_ra)),
  no2_5y_z = seq(-2, 2, by = 0.1),
  pm25_5y_z = 0,
  age_years = mean(analysis_ready$age_years, na.rm = TRUE),
  sex_category = "Female",
  race_category = "White"
)

newdat_no2$pred <- predict(m_int_no2, newdata = newdat_no2, type = "response")

p_int_no2 <- ggplot(newdat_no2, aes(x = no2_5y_z, y = pred, color = traj_cluster_ra)) +
  geom_line(linewidth = 1.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Interaction: NO2 × Trajectory cluster (predicted mortality)",
    x = "NO2 (5y z-score)",
    y = "Predicted probability of death/hospice",
    color = "Cluster"
  ) +
  theme_minimal(base_size = 14)

p_int_no2
save_plot(p_int_no2, "interaction_no2_by_cluster_predprob", w = 10, h = 7)

contrast_pm <- newdat_pm %>%
  filter(pm25_5y_z %in% c(-1, 1)) %>%
  mutate(level = ifelse(pm25_5y_z == -1, "Low (-1 SD)", "High (+1 SD)")) %>%
  select(traj_cluster_ra, level, pred) %>%
  pivot_wider(names_from = level, values_from = pred) %>%
  mutate(
    risk_diff = `High (+1 SD)` - `Low (-1 SD)`,
    risk_ratio = `High (+1 SD)` / `Low (-1 SD)`
  )

p_diff_pm <- ggplot(contrast_pm, aes(x = traj_cluster_ra, y = risk_diff)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    title = "PM2.5 interaction summarized: Δ predicted mortality (+1 SD vs -1 SD)",
    x = "Trajectory cluster",
    y = "Absolute risk difference"
  ) +
  theme_minimal(base_size = 14)

p_diff_pm
save_plot(p_diff_pm, "interaction_pm25_riskdiff_bar", w = 9, h = 6)
save_csv(contrast_pm, "interaction_pm25_contrast_table")

contrast_no2 <- newdat_no2 %>%
  filter(no2_5y_z %in% c(-1, 1)) %>%
  mutate(level = ifelse(no2_5y_z == -1, "Low (-1 SD)", "High (+1 SD)")) %>%
  select(traj_cluster_ra, level, pred) %>%
  pivot_wider(names_from = level, values_from = pred) %>%
  mutate(
    risk_diff = `High (+1 SD)` - `Low (-1 SD)`,
    risk_ratio = `High (+1 SD)` / `Low (-1 SD)`
  )

p_diff_no2 <- ggplot(contrast_no2, aes(x = traj_cluster_ra, y = risk_diff)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    title = "NO2 interaction summarized: Δ predicted mortality (+1 SD vs -1 SD)",
    x = "Trajectory cluster",
    y = "Absolute risk difference"
  ) +
  theme_minimal(base_size = 14)

p_diff_no2
save_plot(p_diff_no2, "interaction_no2_riskdiff_bar", w = 9, h = 6)
save_csv(contrast_no2, "interaction_no2_contrast_table")

heat_pm <- newdat_pm %>%
  mutate(pm_bin = round(pm25_5y_z, 1))

p_heat_pm <- ggplot(heat_pm, aes(x = pm_bin, y = traj_cluster_ra, fill = pred)) +
  geom_tile() +
  scale_fill_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Predicted mortality heatmap: PM2.5 × cluster",
    x = "PM2.5 (5y z-score)",
    y = "Trajectory cluster",
    fill = "Pred. death"
  ) +
  theme_minimal(base_size = 14)

p_heat_pm
save_plot(p_heat_pm, "interaction_pm25_heatmap", w = 10, h = 5.5)

## adjust for confounders 

# Always:
#   
# - age, sex, race/ethnicity
# 
# - year of ICU admission (or spline of date)
# 
# - site (if pooled)
# 
# Add:
#   
# - SES proxy: tract-level deprivation index or median income
# 
# - smoking/lung disease proxy: COPD POA (± asthma/ILD)
# 
# - comorbidity burden: Charlson-like count from POA dx (or just CHF/CKD/DM)
# 
# - cancer severity proxy: metastatic/advanced cancer POA flag


## expand follow up window to 1 week

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(comorbidity)
})

dx_charlson <- hospital_dx %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    diagnosis_code = as.character(diagnosis_code),
    diagnosis_code_format = toupper(coalesce(as.character(diagnosis_code_format), "")),
    poa_present = suppressWarnings(as.integer(poa_present))
  ) %>%
  filter(
    poa_present == 1,
    !is.na(diagnosis_code),
    diagnosis_code != ""
  ) %>%
  mutate(
    code_clean = gsub("[^A-Za-z0-9]", "", diagnosis_code),
    icd_version = case_when(
      str_detect(diagnosis_code_format, "10") ~ "icd10",
      str_detect(diagnosis_code_format, "9")  ~ "icd9",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(icd_version))

charlson_icd10 <- dx_charlson %>%
  filter(icd_version == "icd10") %>%
  select(id = hospitalization_id, code = code_clean) %>%
  comorbidity(
    id = "id",
    code = "code",
    map = "charlson_icd10_quan",
    assign0 = FALSE
  )

charlson_icd10$charlson_score <- score(charlson_icd10, weights = "charlson", assign0 = FALSE)

charlson_all <- charlson_icd10 %>% select(id, charlson_score) %>%
  group_by(id) %>%
  summarize(
    charlson_score = max(charlson_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(hospitalization_id = id)

analysis_ready <- analysis_ready %>%
  mutate(hospitalization_id = as.character(hospitalization_id)) %>%
  left_join(
    charlson_all %>% select(hospitalization_id, charlson_score),
    by = "hospitalization_id"
  )

cluster_charlson <- analysis_ready %>%
  filter(!is.na(traj_cluster_ra)) %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    n = n(),
    charlson_mean = mean(charlson_score, na.rm = TRUE),
    charlson_median = median(charlson_score, na.rm = TRUE),
    charlson_sd = sd(charlson_score, na.rm = TRUE),
    .groups = "drop"
  )

save_csv(cluster_charlson, "cluster_charlson_summary")


 ## lung cancer severity

metastatic_prefixes_icd10 <- c("C77", "C78", "C79", "C80")
metastatic_prefixes_icd9  <- c("196", "197", "198", "199")

pleural_effusion_icd10 <- c("J91", "C785")   # refine as needed
pleural_effusion_icd9  <- c("51181")

cachexia_icd10 <- c("R64")
cachexia_icd9  <- c("7994")

advanced_cancer_flags <- hospital_dx %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    diagnosis_code = norm_code(diagnosis_code),
    diagnosis_code_format = toupper(coalesce(as.character(diagnosis_code_format), "")),
    poa_present = suppressWarnings(as.integer(poa_present)),
    icd10 = str_detect(diagnosis_code_format, "10"),
    icd9  = str_detect(diagnosis_code_format, "9")
  ) %>%
  filter(poa_present == 1) %>%
  group_by(hospitalization_id) %>%
  summarize(
    metastatic_cancer_poa = any(
      (icd10 & code_matches_any_prefix(diagnosis_code, metastatic_prefixes_icd10)) |
        (icd9  & code_matches_any_prefix(diagnosis_code, metastatic_prefixes_icd9)),
      na.rm = TRUE
    ),
    malignant_pleural_effusion_poa = any(
      (icd10 & code_matches_any_prefix(diagnosis_code, pleural_effusion_icd10)) |
        (icd9  & code_matches_any_prefix(diagnosis_code, pleural_effusion_icd9)),
      na.rm = TRUE
    ),
    cachexia_poa = any(
      (icd10 & code_matches_any_prefix(diagnosis_code, cachexia_icd10)) |
        (icd9  & code_matches_any_prefix(diagnosis_code, cachexia_icd9)),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    advanced_cancer_any_poa =
      metastatic_cancer_poa |
      malignant_pleural_effusion_poa |
      cachexia_poa
  )

analysis_ready <- analysis_ready %>%
  left_join(advanced_cancer_flags, by = "hospitalization_id")

analysis_ready %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    metastatic_pct = mean(metastatic_cancer_poa, na.rm = TRUE),
    pleural_effusion_pct = mean(malignant_pleural_effusion_poa, na.rm = TRUE),
    cachexia_pct = mean(cachexia_poa, na.rm = TRUE),
    advanced_any_pct = mean(advanced_cancer_any_poa, na.rm = TRUE),
    .groups = "drop"
  )

# ================================================================================================
# First-24h SOFA for lung cancer ICU trajectory cohort
# Uses existing utils/sofa_calculator.R
# ================================================================================================

# Load SOFA calculator
source(file.path("utils", "sofa_calculator.R"))

# ---------------------------
# 1) ICU admission times for this cohort
# ---------------------------
icu_admit_times <- icu_segments %>%
  semi_join(cohort_lung, by = "hospitalization_id") %>%
  group_by(hospitalization_id) %>%
  summarize(
    icu_admit_time = min(icu_in_time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    icu_admit_time = safe_posix(icu_admit_time)
  )

# ---------------------------
# 2) SOFA cohort input
# ---------------------------
sofa_cohort <- cohort_lung %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id)
  ) %>%
  inner_join(icu_admit_times, by = "hospitalization_id")

# ---------------------------
# 3) Prepare CLIF inputs for calculator
# ---------------------------

# Vitals
vitals_df <- clif_tables[["clif_vitals"]] %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm = safe_posix(recorded_dttm)
  )

# Labs
labs_df <- clif_tables[["clif_labs"]] %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    lab_result_dttm = safe_posix(lab_result_dttm)
  )

# Respiratory / support table
support_df <- rs_raw %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm = safe_posix(recorded_dttm)
  )

# Continuous med admin (for vasoactives in SOFA CV domain)
med_admin_df <- get_min("medication_admin_continuous", c(
  "hospitalization_id",
  "admin_dttm",
  "med_group",
  "med_category",
  "mar_action_group",
  "med_name",
  "med_dose",
  "med_dose_unit"
)) %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    admin_dttm = safe_posix(admin_dttm)
  )

# Patient assessments / scores (for CNS / GCS if calculator uses them)
scores_df <- clif_tables[["clif_patient_assessments"]] %>%
 rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm = safe_posix(recorded_dttm),
    assessment_category = as.character(assessment_category),
    numerical_value = as.numeric(numerical_value)
  ) %>%
  filter(
    !is.na(recorded_dttm) & assessment_category == "gcs_total"
  ) %>%
  select(
    hospitalization_id,
    recorded_dttm,
    assessment_category,
    numerical_value
  )

# ---------------------------
# 4) Safe timestamp helper
# ---------------------------
safe_ts <- function(x) {
  safe_posix(x)
}

# ---------------------------
# 5) Calculate SOFA for first 24h after ICU admit
# ---------------------------
sofa_scores <- calculate_sofa(
  cohort_data  = sofa_cohort,
  vitals_df    = vitals_df,
  labs_df      = labs_df,
  support_df   = support_df,
  med_admin_df = med_admin_df,
  scores_df    = scores_df,
  window_hours = 24,
  safe_ts      = safe_ts
)

# ---------------------------
# 6) Join SOFA into analysis dataset
# ---------------------------
analysis_ready <- analysis_ready %>%
  mutate(hospitalization_id = as.character(hospitalization_id)) %>%
  left_join(
    sofa_scores %>%
      transmute(
        hospitalization_id = as.character(hospitalization_id),
        sofa_total,
        sofa_cv,
        sofa_coag,
        sofa_liver,
        sofa_renal,
        sofa_resp,
        sofa_cns
      ),
    by = "hospitalization_id"
  ) %>%
  mutate(
    across(starts_with("sofa_"), ~ coalesce(.x, 0))
  )

# ---------------------------
# 7) Quick SOFA summary
# ---------------------------
cat("\nSOFA score summary for lung cancer ICU cohort:\n")
print(
  analysis_ready %>%
    select(starts_with("sofa_")) %>%
    summary()
)

# ---------------------------
# 8) SOFA by trajectory cluster
# ---------------------------
cluster_sofa <- analysis_ready %>%
  filter(!is.na(traj_cluster_ra)) %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    n = n(),
    sofa_total_mean   = mean(sofa_total, na.rm = TRUE),
    sofa_total_median = median(sofa_total, na.rm = TRUE),
    sofa_resp_mean    = mean(sofa_resp, na.rm = TRUE),
    sofa_cv_mean      = mean(sofa_cv, na.rm = TRUE),
    sofa_renal_mean   = mean(sofa_renal, na.rm = TRUE),
    sofa_liver_mean   = mean(sofa_liver, na.rm = TRUE),
    sofa_coag_mean    = mean(sofa_coag, na.rm = TRUE),
    sofa_cns_mean     = mean(sofa_cns, na.rm = TRUE),
    .groups = "drop"
  )

save_csv(cluster_sofa, "cluster_sofa_summary")


analysis_ready %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    n = n(),
    sofa_mean = mean(sofa_total, na.rm = TRUE),
    sofa_resp = mean(sofa_resp, na.rm = TRUE),
    sofa_cv = mean(sofa_cv, na.rm = TRUE),
    charlson = mean(charlson_score, na.rm = TRUE),
    advanced_cancer = mean(advanced_cancer_any_poa, na.rm = TRUE),
    mortality = mean(death_or_hospice, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(analysis_ready, aes(traj_cluster_ra, sofa_total)) +
  geom_boxplot() +
  labs(
    x = "Trajectory cluster",
    y = "First 24h SOFA score",
    title = "Baseline ICU severity across respiratory trajectory clusters"
  ) +
  theme_pub()

analysis_ready$admit_year <- year(analysis_ready$admission_dttm)


m_mort_adj <- glm(
  death_or_hospice ~ pm25_5y_z + no2_5y_z + traj_cluster_ra +
    age_years + sex_category + race_category +
    admit_year + charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready,
  family = binomial()
)

summary(m_mort_adj)

m_int_no2_adj <- glm(
  death_or_hospice ~ no2_5y_z * traj_cluster_ra + pm25_5y_z +
    age_years + sex_category + race_category +
    admit_year + charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready,
  family = binomial()
)

summary(m_int_no2_adj)

analysis_ready %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    n = n(),
    sofa_total = mean(sofa_total, na.rm = TRUE),
    sofa_resp = mean(sofa_resp, na.rm = TRUE),
    sofa_cv = mean(sofa_cv, na.rm = TRUE),
    sofa_renal = mean(sofa_renal, na.rm = TRUE),
    sofa_liver = mean(sofa_liver, na.rm = TRUE),
    sofa_coag = mean(sofa_coag, na.rm = TRUE),
    sofa_cns = mean(sofa_cns, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(sofa_cns))

ggplot(analysis_ready, aes(traj_cluster_ra, sofa_cns)) +
  geom_boxplot() +
  labs(
    x = "Trajectory cluster",
    y = "SOFA CNS score",
    title = "CNS organ dysfunction by trajectory cluster"
  ) +
  theme_pub()

gcs_cluster <- scores_df %>%
  left_join(
    analysis_ready %>% select(hospitalization_id, traj_cluster_ra),
    by = "hospitalization_id"
  )

gcs_cluster %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    median_gcs = median(numerical_value, na.rm = TRUE),
    mean_gcs = mean(numerical_value, na.rm = TRUE),
    pct_gcs_lt8 = mean(numerical_value < 8, na.rm = TRUE),
    pct_gcs_lt13 = mean(numerical_value < 13, na.rm = TRUE),
    .groups = "drop"
  )

analysis_ready %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    imv_rate = mean(any_imv_72h, na.rm = TRUE),
    imv_hours = mean(imv_hours_72h, na.rm = TRUE),
    sofa_cns = mean(sofa_cns, na.rm = TRUE),
    .groups = "drop"
  )

neuro_prefix <- c(
  "I61", "I62", "I63", "I64",   # stroke
  "G93",                       # encephalopathy
  "C79.3"                      # brain metastasis
)

neuro_flags <- hospital_dx %>%
  filter(poa_present == 1) %>%
  mutate(code = norm_code(diagnosis_code)) %>%
  filter(code_matches_any_prefix(code, neuro_prefix)) %>%
  distinct(hospitalization_id) %>%
  mutate(neuro_dx = TRUE)

analysis_ready %>%
  left_join(neuro_flags, by = "hospitalization_id") %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    neuro_dx_pct = mean(neuro_dx, na.rm = TRUE),
    .groups = "drop"
  )

cluster_domains <- analysis_ready %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    resp = mean(sofa_resp, na.rm = TRUE),
    cv = mean(sofa_cv, na.rm = TRUE),
    renal = mean(sofa_renal, na.rm = TRUE),
    liver = mean(sofa_liver, na.rm = TRUE),
    coag = mean(sofa_coag, na.rm = TRUE),
    cns = mean(sofa_cns, na.rm = TRUE)
  ) %>%
  pivot_longer(-traj_cluster_ra)

ggplot(cluster_domains, aes(name, value, fill = traj_cluster_ra)) +
  geom_col(position = "dodge") +
  labs(
    x = "SOFA domain",
    y = "Mean score",
    title = "Organ dysfunction signature by trajectory cluster"
  ) +
  theme_pub()


cluster_severity <- analysis_ready %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    mortality = mean(death_or_hospice, na.rm=TRUE),
    sofa = mean(sofa_total, na.rm=TRUE),
    imv_rate = mean(any_imv_72h, na.rm=TRUE),
    icu_los = median(icu_los_hours, na.rm=TRUE),
    .groups="drop"
  )

cluster_severity

cluster_severity <- cluster_severity %>%
  mutate(
    severity_rank = rank(
      mortality + sofa/10 + imv_rate,
      ties.method="first"
    )
  )

library(MASS)

analysis_ready <- analysis_ready %>%
  left_join(
    cluster_severity %>%
      dplyr::select(traj_cluster_ra, severity_rank),
    by="traj_cluster_ra"
  )

analysis_ready$severity_rank <- factor(
  analysis_ready$severity_rank,
  ordered=TRUE
)

m_severity <- polr(
  severity_rank ~ pm25_5y_z + no2_5y_z +
    age_years + sex_category + race_category +
    charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready,
  Hess=TRUE
)

summary(m_severity)

m_traj_adj <- nnet::multinom(
  traj_cluster_ra ~ pm25_5y_z + no2_5y_z +
    age_years + sex_category + race_category +
    charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready
)

summary(m_traj_adj)

m_int_adj <- glm(
  death_or_hospice ~ no2_5y_z * traj_cluster_ra +
    pm25_5y_z +
    age_years + sex_category + race_category +
    charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready,
  family = binomial()
)

summary(m_int_adj)

library(survival)

landmark72 <- analysis_ready %>%
  mutate(
    t0_72 = admission_dttm + dhours(72),
    eligible72 = !is.na(discharge_dttm) & discharge_dttm > t0_72,
    event_after72 = death_or_hospice & discharge_dttm > t0_72,
    time_from72_h = as.numeric(difftime(discharge_dttm, t0_72, units = "hours"))
  ) %>%
  filter(eligible72, !is.na(traj_cluster_ra), !is.na(time_from72_h), time_from72_h >= 0)

cox72 <- coxph(
  Surv(time_from72_h, event_after72) ~
    traj_cluster_ra +
    pm25_5y_z + no2_5y_z +
    age_years + sex_category + race_category +
    charlson_score + sofa_total + advanced_cancer_any_poa,
  data = landmark72
)

summary(cox72)

cox72_int <- coxph(
  Surv(time_from72_h, event_after72) ~ traj_cluster_ra * no2_5y_z + pm25_5y_z +
    age_years + sex_category + race_category + charlson_score + sofa_total + advanced_cancer_any_poa,
  data = landmark72
)

summary(cox72_int)

cor.test(analysis_ready$pm25_5y, analysis_ready$sofa_total)
cor.test(analysis_ready$no2_5y, analysis_ready$sofa_total)

analysis_ready %>%
  group_by(traj_cluster_ra) %>%
  summarize(
    gcs = mean(sofa_cns),
    imv = mean(any_imv_72h),
    brain_met = mean(metastatic_cancer_poa)
  )

# ---------------------------------
# Predicted mortality by cluster x NO2
# ---------------------------------

library(tidyverse)
library(ggplot2)

# Make sure factor levels are stable
analysis_ready <- analysis_ready %>%
  mutate(
    traj_cluster_ra = droplevels(factor(traj_cluster_ra))
  )

# Set reference covariate profile
ref_age <- mean(analysis_ready$age_years, na.rm = TRUE)
ref_charlson <- mean(analysis_ready$charlson_score, na.rm = TRUE)
ref_sofa <- mean(analysis_ready$sofa_total, na.rm = TRUE)

# choose common categories present in your data
ref_sex <- analysis_ready %>%
  count(sex_category, sort = TRUE) %>%
  slice(1) %>%
  pull(sex_category)

ref_race <- analysis_ready %>%
  count(race_category, sort = TRUE) %>%
  slice(1) %>%
  pull(race_category)

ref_adv <- FALSE
if ("advanced_cancer_any_poa" %in% names(analysis_ready)) {
  ref_adv <- FALSE
}

newdat_no2 <- tidyr::expand_grid(
  traj_cluster_ra = levels(analysis_ready$traj_cluster_ra),
  no2_5y_z = seq(-2, 2, by = 0.1)
) %>%
  mutate(
    pm25_5y_z = 0,
    age_years = ref_age,
    sex_category = ref_sex,
    race_category = ref_race,
    charlson_score = ref_charlson,
    sofa_total = ref_sofa,
    advanced_cancer_any_poa = ref_adv
  )

# predicted probability
newdat_no2 <- newdat_no2 %>%
  mutate(
    pred_prob = predict(m_int_adj, newdata = ., type = "response")
  )

p_no2_cluster <- ggplot(newdat_no2,
                        aes(x = no2_5y_z, y = pred_prob, color = traj_cluster_ra)) +
  geom_line(linewidth = 1.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Predicted mortality by trajectory cluster and NO2 exposure",
    x = "NO2 (5-year exposure, z-score)",
    y = "Predicted probability of in-hospital death/hospice",
    color = "Trajectory cluster"
  ) +
  theme_pub()

p_no2_cluster
save_plot(p_no2_cluster, "predicted_mortality_by_cluster_no2", w = 10, h = 7)

# ---------------------------------
# Contrast: NO2 = -1 SD vs +1 SD
# ---------------------------------

contrast_no2 <- newdat_no2 %>%
  filter(no2_5y_z %in% c(-1, 1)) %>%
  mutate(no2_level = ifelse(no2_5y_z == -1, "low", "high")) %>%
  select(traj_cluster_ra, no2_level, pred_prob) %>%
  pivot_wider(names_from = no2_level, values_from = pred_prob) %>%
  mutate(
    risk_diff = high - low,
    risk_ratio = high / low
  )

contrast_no2
save_csv(contrast_no2, "contrast_no2_low_vs_high_by_cluster")

p_no2_diff <- ggplot(contrast_no2, aes(x = traj_cluster_ra, y = risk_diff, fill = traj_cluster_ra)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    title = "Absolute mortality increase from low to high NO2 by trajectory cluster",
    x = "Trajectory cluster",
    y = "Risk difference: high NO2 - low NO2",
    fill = "Trajectory cluster"
  ) +
  theme_pub() +
  theme(legend.position = "none")

p_no2_diff
save_plot(p_no2_diff, "risk_difference_no2_by_cluster", w = 9, h = 6)


m_total_no2 <- glm(
  death_or_hospice ~ no2_5y_z + pm25_5y_z +
    age_years + sex_category + race_category +
    charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready,
  family = binomial()
)

summary(m_total_no2)

m_direct_no2 <- glm(
  death_or_hospice ~ no2_5y_z + pm25_5y_z + severity_rank +
    age_years + sex_category + race_category +
    charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready,
  family = binomial()
)

summary(m_direct_no2)

beta_total_no2  <- coef(m_total_no2)[["no2_5y_z"]]
beta_direct_no2 <- coef(m_direct_no2)[["no2_5y_z"]]

mediation_decomp <- tibble(
  exposure = "NO2_5y_z",
  beta_total = beta_total_no2,
  beta_direct = beta_direct_no2,
  beta_indirect_approx = beta_total_no2 - beta_direct_no2,
  pct_mediated_approx = (beta_total_no2 - beta_direct_no2) / beta_total_no2
)

mediation_decomp
save_csv(mediation_decomp, "mediation_style_decomposition_no2")

m_total_no2_cluster <- glm(
  death_or_hospice ~ no2_5y_z + pm25_5y_z +
    age_years + sex_category + race_category +
    charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready,
  family = binomial()
)

m_direct_no2_cluster <- glm(
  death_or_hospice ~ no2_5y_z + pm25_5y_z + traj_cluster_ra +
    age_years + sex_category + race_category +
    charlson_score + sofa_total + advanced_cancer_any_poa,
  data = analysis_ready,
  family = binomial()
)

beta_total_no2_cluster  <- coef(m_total_no2_cluster)[["no2_5y_z"]]
beta_direct_no2_cluster <- coef(m_direct_no2_cluster)[["no2_5y_z"]]

mediation_decomp_cluster <- tibble(
  exposure = "NO2_5y_z",
  beta_total = beta_total_no2_cluster,
  beta_direct = beta_direct_no2_cluster,
  beta_indirect_approx = beta_total_no2_cluster - beta_direct_no2_cluster,
  pct_mediated_approx = (beta_total_no2_cluster - beta_direct_no2_cluster) / beta_total_no2_cluster
)

mediation_decomp_cluster
save_csv(mediation_decomp_cluster, "mediation_style_decomposition_no2_cluster")

analysis_ready <- analysis_ready %>%
  mutate(
    no2_q = ntile(no2_5y, 4)
  )

p_sev_no2 <- ggplot(analysis_ready, aes(x = factor(no2_q), fill = severity_rank)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Trajectory severity distribution across NO2 quartiles",
    x = "NO2 quartile",
    y = "Proportion",
    fill = "Severity rank"
  ) +
  theme_pub()

p_sev_no2
save_plot(p_sev_no2, "severity_distribution_by_no2_quartile", w = 9, h = 6)















