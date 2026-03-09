





suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(TraMineR)
  library(cluster)
})

# --------------------------- Trajectory: hourly device state (0–72h) ---------------------------


H_MAX <- 72L

keys <- cohort_lung %>%
  transmute(hospitalization_id = as.character(hospitalization_id),
            t0 = as.POSIXct(t0, tz="UTC"))

rs_for_traj <- rs_raw %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm      = as.POSIXct(recorded_dttm, tz="UTC"),
    device_category    = as.character(device_category)
  ) %>%
  inner_join(keys, by="hospitalization_id") %>%
  mutate(
    dt_h = as.numeric(difftime(recorded_dttm, t0, units="hours")),
    h    = floor(dt_h),
    dev  = device_state_clif(device_category),
    rank = device_rank_clif(device_category)
  ) %>%
  filter(!is.na(h), h >= 0, h <= H_MAX, !is.na(rank))

# Dominant device per hour = max rank; if ties, take first after arranging by recorded_dttm
hourly_state <- rs_for_traj %>%
  arrange(hospitalization_id, h, desc(rank), recorded_dttm) %>%
  group_by(hospitalization_id, h) %>%
  slice(1) %>%
  ungroup() %>%
  select(hospitalization_id, h, state = dev)

# Complete 0..H_MAX and fill missing with LOCF/NOCB
panel <- hourly_state %>%
  group_by(hospitalization_id) %>%
  complete(h = 0:H_MAX) %>%
  arrange(hospitalization_id, h) %>%
  tidyr::fill(state, .direction = "down") %>%   # LOCF
  tidyr::fill(state, .direction = "up") %>%     # NOCB
  mutate(state = coalesce(state, "UNK")) %>%
  ungroup()

# Wide
seq_wide <- panel %>%
  mutate(h_col = paste0("H", h)) %>%
  select(hospitalization_id, h_col, state) %>%
  pivot_wider(names_from = h_col, values_from = state)

# Alphabet in clinical order
alphabet <- c("ROOM AIR","NASAL CANNULA","FACE MASK","TRACH COLLAR","HIGH FLOW NC","CPAP","NIPPV","IMV","OTHER","UNK")

seq_cols <- seq_wide %>% select(matches("^H\\d+$"))

seq_obj <- TraMineR::seqdef(
  seq_cols,
  alphabet = alphabet,
  states   = alphabet,
  labels   = alphabet,
  xtstep   = 6
)

# Distance + clustering
dist_om <- TraMineR::seqdist(seq_obj, method = "OM", indel = 1, sm = "TRATE")
hc <- hclust(as.dist(dist_om), method = "ward.D2")

K <- 5
cluster <- cutree(hc, k = K)

traj_assign <- tibble(
  hospitalization_id = seq_wide$hospitalization_id,
  traj_cluster = as.integer(cluster)
)


# Plots (base graphics)
png(file.path(out_dir, make_name("traj_seqdplot_by_cluster", "png")), width=1400, height=1400, res=160)
TraMineR::seqdplot(seq_obj, group = cluster, border = NA, main = "CLIF respiratory support distribution by cluster (0–72h)")
dev.off()

png(file.path(out_dir, make_name("traj_seqrplot_by_cluster", "png")), width=1500, height=1700, res=160)
TraMineR::seqrplot(seq_obj, group = cluster, method = "OM", indel = 1, sm = "TRATE",
                   main = "Representative CLIF respiratory support sequences by cluster")
dev.off()

# Summaries
traj_summary <- analysis_ready %>%
  left_join(traj_assign, by = "hospitalization_id") %>%
  group_by(traj_cluster) %>%
  summarize(
    n = n(),
    death_or_hospice = mean(death_or_hospice, na.rm=TRUE),
    icu_los_h_median = median(icu_los_hours, na.rm=TRUE),
    any_imv_72h = mean(any_imv_72h, na.rm=TRUE),
    imv_hours_72h_median = median(imv_hours_72h, na.rm=TRUE),
    pm25_3y_median = median(pm25_3y, na.rm=TRUE),
    no2_3y_median  = median(no2_3y,  na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(desc(n))

save_csv(traj_summary, "traj_cluster_summary")

med_cont <- get_min("medication_admin_continuous", c(
  "hospitalization_id",
  "admin_dttm",
  "med_group",
  "mar_action_group"
)) %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    admin_dttm         = safe_posix(admin_dttm),
    med_group          = tolower(as.character(med_group)),
    mar_action_group   = tolower(as.character(mar_action_group))
  )

vaso_admin <- med_cont %>%
  filter(
    med_group == "vasoactives",
    mar_action_group == "administered",
    !is.na(admin_dttm)
  )

keys <- cohort_lung %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    t0 = as.POSIXct(t0, tz="UTC")
  )

H_MAX <- 72L

vaso_hourly <- vaso_admin %>%
  inner_join(keys, by="hospitalization_id") %>%
  mutate(
    dt_h = as.numeric(difftime(admin_dttm, t0, units="hours")),
    h    = floor(dt_h)
  ) %>%
  filter(!is.na(h), h >= 0, h <= H_MAX) %>%
  distinct(hospitalization_id, h) %>%
  mutate(vaso_h = 1)

vitals <- get_min("vitals", c(
  "hospitalization_id","recorded_dttm","vital_category","vital_value"
)) %>%
  rename_with(tolower) %>%
  mutate(
    recorded_dttm   = safe_posix(recorded_dttm),
    vital_category  = tolower(as.character(vital_category)),
    vital_value     = as_num(vital_value)
  )

H_MAX <- 72L

keys <- cohort_lung %>%
  transmute(hospitalization_id = as.character(hospitalization_id),
            t0 = as.POSIXct(t0, tz="UTC"))

vitals_hourly <- vitals %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm,
    vital_category,
    vital_value
  ) %>%
  inner_join(keys, by="hospitalization_id") %>%
  mutate(
    dt_h = as.numeric(difftime(recorded_dttm, t0, units="hours")),
    h = floor(dt_h)
  ) %>%
  filter(!is.na(h), h >= 0, h <= H_MAX) %>%
  # keep only categories we care about for the trajectory
  filter(vital_category %in% c("spo2","respiratory_rate","map","heart_rate")) %>%
  group_by(hospitalization_id, h, vital_category) %>%
  summarize(v = suppressWarnings(median(vital_value, na.rm=TRUE)), .groups="drop") %>%
  mutate(v = ifelse(is.infinite(v), NA, v)) %>%
  tidyr::pivot_wider(names_from = vital_category, values_from = v) %>%
  rename(
    spo2_h = spo2,
    rr_h   = respiratory_rate,
    map_h  = map,
    hr_h   = heart_rate
  )

# assumes device_state_clif() and device_rank_clif() are defined
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
    dev  = device_state_clif(device_category),
    rank = device_rank_clif(device_category)
  ) %>%
  filter(!is.na(h), h >= 0, h <= H_MAX, !is.na(rank))

hourly_state <- rs_for_traj %>%
  arrange(hospitalization_id, h, desc(rank), recorded_dttm) %>%
  group_by(hospitalization_id, h) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(hospitalization_id, h, rs_state = dev)

rs_panel <- hourly_state %>%
  group_by(hospitalization_id) %>%
  complete(h = 0:H_MAX) %>%
  arrange(hospitalization_id, h) %>%
  tidyr::fill(rs_state, .direction = "down") %>%
  tidyr::fill(rs_state, .direction = "up") %>%
  mutate(rs_state = coalesce(rs_state, "UNK")) %>%
  ungroup()


phys_panel <- rs_panel %>%
  left_join(vaso_hourly, by=c("hospitalization_id","h")) %>%
  mutate(
    vaso_h = coalesce(vaso_h, 0)
  )

phys_panel %>%
  group_by(hospitalization_id) %>%
  summarize(
    any_vaso = max(vaso_h),
    total_hours_vaso = sum(vaso_h)
  ) %>%
  summarize(
    pct_any_vaso = mean(any_vaso),
    mean_vaso_hours = mean(total_hours_vaso)
  )


### ARF definition

vitals <- clif_tables[["clif_vitals"]]
labs <- clif_tables[["clif_labs"]]

H_MAX <- 72L
ROOM_AIR_FIO2 <- 0.21
JOIN_NEAR_H   <- 1
HYPERPAIR_H   <- 2

# --- keys: hosp + t0 ---
keys <- cohort_lung %>%
  transmute(
    hospitalization_id = as.character(hospitalization_id),
    t0 = as.POSIXct(t0, tz="UTC")
  )

# --- pull signals in [t0, t0+72h] (trajectory window) ---
win72 <- keys %>%
  transmute(
    hospitalization_id,
    win_start = t0,
    win_end   = t0 + dhours(H_MAX)
  )

# SpO2 from CLIF vitals long
spo2_win <- vitals %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm = safe_posix(recorded_dttm),
    vital_category = tolower(as.character(vital_category)),
    vital_value = as_num(vital_value)
  ) %>%
  filter(vital_category == "spo2") %>%
  inner_join(win72, by="hospitalization_id") %>%
  filter(recorded_dttm >= win_start, recorded_dttm <= win_end) %>%
  transmute(hospitalization_id, spo2_time = recorded_dttm, spo2 = vital_value)

# ABG labs long
labs_win <- labs %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    lab_result_dttm = safe_posix(lab_result_dttm),
    lab_category = tolower(as.character(lab_category)),
    lab_value_numeric = as_num(lab_value_numeric)
  ) %>%
  filter(lab_category %in% c("po2_arterial","pco2_arterial","ph_arterial")) %>%
  inner_join(win72, by="hospitalization_id") %>%
  filter(lab_result_dttm >= win_start, lab_result_dttm <= win_end) %>%
  transmute(hospitalization_id, lab_time = lab_result_dttm, lab_category, val = lab_value_numeric)

# FiO2 from respiratory support (must have fio2_set)
fio2_win <- rs_raw %>%
  rename_with(tolower) %>%
  mutate(
    hospitalization_id = as.character(hospitalization_id),
    recorded_dttm = safe_posix(recorded_dttm),
    fio2_set = as_num(fio2_set)
  ) %>%
  inner_join(win72, by="hospitalization_id") %>%
  filter(recorded_dttm >= win_start, recorded_dttm <= win_end) %>%
  transmute(hospitalization_id, fio2_time = recorded_dttm, fio2_set)

# --- data.table nearest pairing helpers ---
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

# PaO2 ↔ FiO2 nearest within JOIN_NEAR_H for P/F and RA PaO2 rule
po2dt <- labs_win[lab_category == "po2_arterial", .(hospitalization_id, po2_time = lab_time, po2 = val)]
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
pco2dt <- labs_win[lab_category == "pco2_arterial", .(hospitalization_id, pco2_time = lab_time, pco2 = val)]
phdt   <- labs_win[lab_category == "ph_arterial",   .(hospitalization_id, ph_time   = lab_time, ph   = val)]
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

# --- Convert event timestamps to hour bins relative to t0 ---
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
    # hypox components
    hypox_roomair_spo2 = coalesce(hypox_roomair_spo2, FALSE),
    hypox_roomair_po2  = coalesce(hypox_roomair_po2,  FALSE),
    hypox_pf           = coalesce(hypox_pf,           FALSE),
    
    # hyper component
    hypercapnic        = coalesce(hypercapnic,        FALSE),
    
    # collapse to two core dimensions
    hypox_h = hypox_roomair_spo2 | hypox_roomair_po2 | hypox_pf,
    hyper_h = hypercapnic
  ) %>%
  select(hospitalization_id, h, hypox_h, hyper_h,
         hypox_roomair_spo2, hypox_roomair_po2, hypox_pf, hypercapnic)


carry_forward <- function(x, carry = 6L) {
  out <- rep(FALSE, length(x))
  idx <- which(x)
  for (i in idx) out[i:min(length(x), i + carry)] <- TRUE
  out
}

arf_hourly <- arf_hourly %>%
  group_by(hospitalization_id) %>%
  complete(
    h = 0:H_MAX,
    fill = list(
      hypox_h = FALSE, hyper_h = FALSE,
      hypox_roomair_spo2 = FALSE, hypox_roomair_po2 = FALSE, hypox_pf = FALSE,
      hypercapnic = FALSE
    )
  ) %>%
  arrange(h) %>%
  mutate(
    hypox_h = carry_forward(hypox_h, carry = 6L),
    hyper_h = carry_forward(hyper_h, carry = 6L)
  ) %>%
  ungroup()

arf_hourly <- arf_hourly %>%
  mutate(
    arf_subtype_h = case_when(
      hypox_h & hyper_h ~ "ARF_MIXED",
      hypox_h & !hyper_h ~ "ARF_HYPOX",
      !hypox_h & hyper_h ~ "ARF_HYPER",
      TRUE ~ "NO_ARF"
    )
  )

panel_mof <- rs_panel %>%
  left_join(arf_hourly %>% select(hospitalization_id, h, arf_subtype_h),
            by = c("hospitalization_id","h")) %>%
  left_join(vaso_hourly, by = c("hospitalization_id","h")) %>%
  mutate(
    arf_subtype_h = coalesce(arf_subtype_h, "NO_ARF"),
    vaso_h = coalesce(vaso_h, 0L),
    state = paste(
      rs_state,
      arf_subtype_h,
      ifelse(vaso_h == 1L, "VASO", "NoVASO"),
      sep = "|"
    )
  )

















analysis_ready <- analysis_ready %>%
  left_join(traj_assign, by = "hospitalization_id")


library(nnet)

analysis_ready <- analysis_ready %>%
  mutate(
    traj_cluster = factor(traj_cluster),
    pm25_5y_z = scale(pm25_5y),
    no2_5y_z  = scale(no2_5y)
  )

m_traj <- multinom(
  traj_cluster ~ pm25_5y_z + no2_5y_z +
    age_years + sex_category + race_category,
  data = analysis_ready
)

summary(m_traj)


m_mort <- glm(
  death_or_hospice ~ pm25_5y_z + no2_5y_z +
    traj_cluster +
    age_years + sex_category + race_category,
  data = analysis_ready,
  family = binomial()
)

summary(m_mort)
exp(cbind(OR = coef(m_mort), confint(m_mort)))

m_int <- glm(
  death_or_hospice ~ pm25_5y_z * traj_cluster +
    no2_5y_z +
    age_years + sex_category + race_category,
  data = analysis_ready,
  family = binomial()
)

summary(m_int)

m_int2 <- glm(
  death_or_hospice ~ no2_5y_z * traj_cluster +
    pm25_5y_z +
    age_years + sex_category + race_category,
  data = analysis_ready,
  family = binomial()
)

summary(m_int2)

m_los <- glm(
  icu_los_hours ~ pm25_5y_z + no2_5y_z +
    traj_cluster +
    age_years + sex_category + race_category,
  data = analysis_ready,
  family = quasipoisson(link = "log")
)

summary(m_los)


























# phys_panel <- rs_panel %>%
#   left_join(vitals_hourly, by=c("hospitalization_id","h")) %>%
#   group_by(hospitalization_id) %>%
#   complete(h = 0:H_MAX) %>%
#   arrange(hospitalization_id, h) %>%
#   # keep RS fully populated
#   tidyr::fill(rs_state, .direction = "downup") %>%
#   ungroup()
# 
# spo2_bin <- function(x) {
#   case_when(
#     is.na(x)      ~ "SPO2_UNK",
#     x < 88        ~ "SPO2_VLOW",
#     x >= 88 & x < 92 ~ "SPO2_LOW",
#     x >= 92 & x < 96 ~ "SPO2_OK",
#     x >= 96       ~ "SPO2_HI"
#   )
# }
# 
# rr_bin <- function(x) {
#   case_when(
#     is.na(x)         ~ "RR_UNK",
#     x < 12           ~ "RR_LOW",
#     x >= 12 & x < 20 ~ "RR_NL",
#     x >= 20 & x < 30 ~ "RR_HI",
#     x >= 30          ~ "RR_VHI"
#   )
# }
# 
# map_bin <- function(x) {
#   case_when(
#     is.na(x)          ~ "MAP_UNK",
#     x < 65            ~ "MAP_LOW",
#     x >= 65 & x < 80  ~ "MAP_OK",
#     x >= 80           ~ "MAP_HI"
#   )
# }
# 
# phys_panel <- phys_panel %>%
#   mutate(
#     spo2_cat = spo2_bin(spo2_h),
#     rr_cat   = rr_bin(rr_h),
#     map_cat  = map_bin(map_h)   # optional: keep if you want hemodynamics in the state
#   )
# 
# phys_panel <- phys_panel %>%
#   mutate(
#     comp_state = paste(rs_state, spo2_cat, rr_cat, sep="|")
#     # if adding MAP: paste(rs_state, spo2_cat, rr_cat, map_cat, sep="|")
#   )
# 
# seq_wide_phys <- phys_panel %>%
#   mutate(h_col = paste0("H", h)) %>%
#   select(hospitalization_id, h_col, comp_state) %>%
#   pivot_wider(names_from = h_col, values_from = comp_state)
# 
# seq_cols_phys <- seq_wide_phys %>% select(matches("^H\\d+$"))
# 
# alphabet_phys <- sort(unique(unlist(seq_cols_phys)))
# 
# seq_obj_phys <- TraMineR::seqdef(
#   seq_cols_phys,
#   alphabet = alphabet_phys,
#   states   = alphabet_phys,
#   labels   = alphabet_phys,
#   xtstep   = 6
# )
# 
# dist_om_phys <- TraMineR::seqdist(seq_obj_phys, method="OM", indel=1, sm="TRATE")
# hc_phys <- hclust(as.dist(dist_om_phys), method="ward.D2")
# 
# K <- 6
# cluster_phys <- cutree(hc_phys, k = K)
# 
# traj_assign_phys <- tibble(
#   hospitalization_id = seq_wide_phys$hospitalization_id,
#   traj_cluster_phys  = as.integer(cluster_phys)
# )
# 
# save_csv(traj_assign_phys, "traj_cluster_assignments_phys")
# 
# marg_rs <- phys_panel %>%
#   left_join(traj_assign_phys, by="hospitalization_id") %>%
#   group_by(traj_cluster_phys, h, rs_state) %>%
#   summarize(p = n() / sum(n()), .groups="drop")
# 
# marg_spo2 <- phys_panel %>%
#   left_join(traj_assign_phys, by="hospitalization_id") %>%
#   group_by(traj_cluster_phys, h, spo2_cat) %>%
#   summarize(p = n() / sum(n()), .groups="drop")
# 
# marg_rr <- phys_panel %>%
#   left_join(traj_assign_phys, by="hospitalization_id") %>%
#   group_by(traj_cluster_phys, h, rr_cat) %>%
#   summarize(p = n() / sum(n()), .groups="drop")
# 
# save_csv(marg_rs,   "traj_phys_marginal_rs_by_cluster_hour")
# save_csv(marg_spo2, "traj_phys_marginal_spo2_by_cluster_hour")
# save_csv(marg_rr,   "traj_phys_marginal_rr_by_cluster_hour")
# 
# p1 <- ggplot(marg_rs, aes(x=h, y=p, fill=rs_state)) +
#   geom_area() +
#   facet_wrap(~traj_cluster_phys, ncol=2) +
#   labs(title="Respiratory support state distribution by cluster (0–72h)", x="Hour", y="Proportion") +
#   theme_pub()
# 
# p1
# 
# p2 <- ggplot(marg_spo2, aes(x=h, y=p, fill=spo2_cat)) +
#   geom_area() +
#   facet_wrap(~traj_cluster_phys, ncol=2) +
#   labs(title="SpO2 state distribution by cluster (0–72h)", x="Hour", y="Proportion") +
#   theme_pub()
# 
# p2
# 
# p3 <- ggplot(marg_rr, aes(x=h, y=p, fill=rr_cat)) +
#   geom_area() +
#   facet_wrap(~traj_cluster_phys, ncol=2) +
#   labs(title="Resp rate state distribution by cluster (0–72h)", x="Hour", y="Proportion") +
#   theme_pub()
# 
# p3
# 
# save_plot(p1, "traj_phys_rs_area", w=11, h=7)
# 

