# ==========================
# County-level exposome linkage (NO2, PM2.5)
# Inputs:
#   - cohort_lung (from your cohort builder)
#   - exposome files: no2_county_year, pm25_county_year in repo/exposome/
# Output:
#   - cohort_lung_expo (cohort_lung + linked chronic exposures)
# ==========================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(lubridate)
  library(arrow)
})

# ---- User parameters ----
EXPO_DIR <- file.path(repo, "exposome")   # adjust if different
NO2_FILE <- file.path(EXPO_DIR, "no2_county_year")
PM25_FILE <- file.path(EXPO_DIR, "pm25_county_year")

# Exposure window: mean over prior N years, excluding index year
EXPO_LAG_YEARS <- 5   # 3 or 5 are common; keep 5 as default
EXCLUDE_INDEX_YEAR <- TRUE

# If you want to anchor exposure to ICU-in year, use t0; otherwise admission_dttm
INDEX_TIME_VAR <- "t0"   # "t0" recommended given your pipeline

# ---- Helpers ----
read_any_local <- function(path_no_ext) {
  # Accepts either exact file with extension or prefix without extension.
  # Tries: .parquet, .csv, .fst
  candidates <- c(
    path_no_ext,
    paste0(path_no_ext, ".parquet"),
    paste0(path_no_ext, ".csv"),
    paste0(path_no_ext, ".fst")
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) stop("No exposome file found for prefix: ", path_no_ext)
  
  path <- candidates[1]
  ext <- tolower(tools::file_ext(path))
  message("Reading exposome: ", path)
  switch(ext,
         "parquet" = arrow::read_parquet(path),
         "csv"     = readr::read_csv(path, show_col_types = FALSE),
         "fst"     = fst::read_fst(path, as.data.table = FALSE),
         stop("Unsupported exposome extension: ", ext))
}

pad5 <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "[^0-9]", "")
  ifelse(nchar(x) == 0, NA_character_, str_pad(x, 5, pad = "0"))
}

# ---- Load exposome tables ----
no2 <- read_any_local(NO2_FILE) %>% rename_with(tolower)
pm25 <- read_any_local(PM25_FILE) %>% rename_with(tolower)

# ---- Standardize expected columns ----
# We try to infer the value column if not explicitly named no2/pm25.
infer_value_col <- function(df, prefer) {
  # prefer: c("no2","pm25")
  nms <- names(df)
  # common patterns
  cand <- nms[str_detect(nms, paste0("^", prefer, "($|_)"))]
  if (length(cand) > 0) return(cand[1])
  # else take first numeric column besides county/year
  num_cand <- nms[sapply(df, is.numeric)]
  num_cand <- setdiff(num_cand, c("year"))
  if (length(num_cand) == 0) stop("Cannot infer value column for ", prefer, " from: ", paste(nms, collapse=", "))
  num_cand[1]
}

# County/year column detection
infer_county_col <- function(df) {
  nms <- names(df)
  cand <- nms[str_detect(nms, "county") & str_detect(nms, "fips|code|geoid")]
  if (length(cand) > 0) return(cand[1])
  # fallback: exact matches
  cand2 <- intersect(nms, c("county_fips","county_code","fips","geoid"))
  if (length(cand2) > 0) return(cand2[1])
  stop("Cannot infer county column from: ", paste(nms, collapse=", "))
}

county_col_no2 <- infer_county_col(no2)
county_col_pm25 <- infer_county_col(pm25)

val_col_no2 <- infer_value_col(no2, "no2")
val_col_pm25 <- infer_value_col(pm25, "pm25")

no2_std <- no2 %>%
  transmute(
    county_fips = pad5(.data[[county_col_no2]]),
    year = as.integer(year),
    no2_value = suppressWarnings(as.numeric(.data[[val_col_no2]]))
  ) %>%
  filter(!is.na(county_fips), !is.na(year))

pm25_std <- pm25 %>%
  transmute(
    county_fips = pad5(.data[[county_col_pm25]]),
    year = as.integer(year),
    pm25_value = suppressWarnings(as.numeric(.data[[val_col_pm25]]))
  ) %>%
  filter(!is.na(county_fips), !is.na(year))

expo_county_year <- no2_std %>%
  full_join(pm25_std, by = c("county_fips","year"))

# ---- Build cohort index year ----
stopifnot(INDEX_TIME_VAR %in% names(cohort_lung))

cohort_geo <- cohort_lung %>%
  mutate(
    county_fips = pad5(county_code),
    index_time = as.POSIXct(.data[[INDEX_TIME_VAR]], tz = "UTC"),
    index_year = lubridate::year(index_time)
  ) %>%
  select(hospitalization_id, county_fips, index_year) %>%
  filter(!is.na(county_fips), !is.na(index_year))

# ---- Compute lagged multi-year mean prior to index year ----
cohort_expo_long <- cohort_geo %>%
  left_join(expo_county_year, by = "county_fips") %>%
  mutate(
    year_start = index_year - EXPO_LAG_YEARS,
    year_end = if (EXCLUDE_INDEX_YEAR) index_year - 1 else index_year
  ) %>%
  filter(year >= year_start, year <= year_end)

expo_agg <- cohort_expo_long %>%
  group_by(hospitalization_id) %>%
  summarise(
    no2_mean_lag = mean(no2_value, na.rm = TRUE),
    pm25_mean_lag = mean(pm25_value, na.rm = TRUE),
    n_years_no2 = sum(!is.na(no2_value)),
    n_years_pm25 = sum(!is.na(pm25_value)),
    .groups = "drop"
  )

# ---- Merge back onto cohort ----
cohort_lung_expo <- cohort_lung %>%
  left_join(expo_agg, by = "hospitalization_id")

# ---- Sanity checks ----
message("Exposome linkage summary:")
print(cohort_lung_expo %>%
        summarise(
          n = n(),
          missing_county = sum(is.na(county_code)),
          missing_expo_no2 = sum(is.na(no2_mean_lag)),
          missing_expo_pm25 = sum(is.na(pm25_mean_lag)),
          mean_no2 = mean(no2_mean_lag, na.rm = TRUE),
          mean_pm25 = mean(pm25_mean_lag, na.rm = TRUE)
        ))

# Optional: drop if too sparse
# cohort_lung_expo <- cohort_lung_expo %>% filter(n_years_no2 >= 3, n_years_pm25 >= 3)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(stringr)
  library(ggplot2)
  library(nnet)        # multinom
  library(survival)    # coxph
  library(broom)
})

# ----------------------------
# Assumes you already have:
#   - cohort_lung_expo  (1 row/hosp; includes proc_time, icu_in_time/t0, county_code, exposures)
#   - clif_tables list with clif_respiratory_support, clif_labs, clif_vitals
# ----------------------------

# ---- Helpers ----
safe_posix <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  if (is.numeric(x)) return(as.POSIXct(x, origin = "1970-01-01", tz = "UTC"))
  suppressWarnings(as.POSIXct(x, tz = "UTC"))
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# Map device_category to ordered respiratory support states
map_state <- function(device_category) {
  d <- as.character(device_category)
  case_when(
    d == "IMV" ~ "IMV",
    d %in% c("NIPPV","CPAP") ~ "NIV",
    d == "High Flow NC" ~ "HFNC",
    d %in% c("NC","Nasal Cannula","Simple Mask","Non-Rebreather","Face Mask","Venturi") ~ "LowO2",
    d %in% c("Room Air","RA") ~ "RA",
    is.na(d) | d == "" ~ "None",
    TRUE ~ "Other"
  )
}
state_levels <- c("None","RA","LowO2","HFNC","NIV","IMV","Other")
state_sev <- setNames(seq_along(state_levels), state_levels)  # higher = “more severe” if ordered

# Identify candidate lab/vital names robustly (sites vary)
pick_first_match <- function(x, patterns) {
  # x: character vector of names; patterns: regex vector
  out <- x[Reduce(`|`, lapply(patterns, function(p) str_detect(tolower(x), p)))]
  if (length(out) == 0) NA_character_ else out[1]
}

# Convert FiO2 to fraction (0-1)
fio2_to_frac <- function(fio2_raw) {
  f <- as_num(fio2_raw)
  case_when(
    is.na(f) ~ NA_real_,
    f > 1.5 ~ f/100,   # assume percent
    TRUE ~ f           # already fraction
  )
}

# SpO2 to proportion (0-1)
spo2_to_prop <- function(spo2_raw) {
  s <- as_num(spo2_raw)
  case_when(
    is.na(s) ~ NA_real_,
    s > 1.5 ~ s/100,
    TRUE ~ s
  )
}

# PaO2: assume mmHg numeric
pao2_to_mmhg <- function(pao2_raw) as_num(pao2_raw)

# ---------------------------------------
# 1) Pull CLIF tables and harmonize times
# ---------------------------------------
rs <- clif_tables[["clif_respiratory_support"]] %>%
  rename_with(tolower) %>%
  mutate(recorded_dttm = safe_posix(recorded_dttm))

labs <- clif_tables[["clif_labs"]] %>%
  rename_with(tolower) %>%
  mutate(
    lab_collect_dttm = safe_posix(lab_collect_dttm),
    lab_result_dttm  = safe_posix(lab_result_dttm),
    lab_value_num = dplyr::coalesce(
      suppressWarnings(as.numeric(lab_value_numeric)),
      suppressWarnings(as.numeric(as.character(lab_value)))
    )
  )


vitals <- clif_tables[["clif_vitals"]] %>%
  rename_with(tolower) %>%
  mutate(recorded_dttm = safe_posix(recorded_dttm))

# Use t0 = ICU in time (you set it that way)
t0_tbl <- cohort_lung_expo %>%
  transmute(
    hospitalization_id,
    proc_time = safe_posix(proc_time),
    t0 = safe_posix(icu_in_time),
    discharge_dttm = safe_posix(discharge_dttm)
  )

# ---------------------------------------
# 2) Respiratory trajectories (0–72h, 6h)
# ---------------------------------------
rs_panel <- rs %>%
  inner_join(t0_tbl %>% select(hospitalization_id, t0), by="hospitalization_id") %>%
  mutate(
    dt_h = as.numeric(difftime(recorded_dttm, t0, units="hours")),
    bin = floor(dt_h / 6) * 6
  ) %>%
  filter(dt_h >= 0, dt_h <= 72) %>%
  mutate(
    state = factor(map_state(device_category), levels = state_levels),
    sev = as.integer(state)  # respects factor order above
  ) %>%
  group_by(hospitalization_id, bin) %>%
  summarise(
    # choose the most severe state within the bin
    state = state[which.max(sev)][1],
    .groups="drop"
  )

# Make complete sequence (fill missing bins as "None")
all_bins <- seq(0, 72, by=6)
rs_seq <- rs_panel %>%
  tidyr::complete(hospitalization_id, bin = all_bins, fill = list(state = factor("None", levels=state_levels))) %>%
  arrange(hospitalization_id, bin)

# Wide matrix
rs_wide <- rs_seq %>%
  mutate(bin_chr = paste0("h", sprintf("%02d", bin))) %>%
  select(hospitalization_id, bin_chr, state) %>%
  pivot_wider(names_from = bin_chr, values_from = state)

# Simple, interpretable classing (high value, low fragility):
# Class 1: early IMV (any IMV in first 24h)
# Class 2: escalation (max severity increases by >=2 levels over 72h)
# Class 3: stable low (never above LowO2)
# Class 4: persistent advanced noninvasive (HFNC/NIV but no IMV)
# (You can replace this later with TraMineR clustering; start here.)

rs_class <- rs_seq %>%
  mutate(sev = as.integer(state)) %>%
  group_by(hospitalization_id) %>%
  summarise(
    max_sev_24 = max(sev[bin <= 24], na.rm=TRUE),
    max_sev_72 = max(sev, na.rm=TRUE),
    min_sev_72 = min(sev, na.rm=TRUE),
    delta_sev  = max_sev_72 - min_sev_72,
    any_imv_24 = any(state == "IMV" & bin <= 24),
    any_imv_72 = any(state == "IMV"),
    any_hfnc   = any(state == "HFNC"),
    any_niv    = any(state == "NIV"),
    class = case_when(
      any_imv_24 ~ "Early IMV",
      !any_imv_72 & (any_hfnc | any_niv) & max_sev_72 >= state_sev[["HFNC"]] ~ "Persistent HFNC/NIV",
      delta_sev >= 2 ~ "Escalation",
      max_sev_72 <= state_sev[["LowO2"]] ~ "Stable low support",
      TRUE ~ "Other"
    ),
    .groups="drop"
  ) %>%
  mutate(class = factor(class, levels = c("Stable low support","Persistent HFNC/NIV","Escalation","Early IMV","Other")))

# ---------------------------------------
# 3) Oxygenation severity (PF and SF)
# ---------------------------------------
# Identify PaO2 lab name and creatinine lab name
lab_names <- unique(na.omit(labs$lab_name))
pao2_name <- pick_first_match(lab_names, c("pao2", "pa o2", "arterial.*o2", "^po2$"))
creat_name <- pick_first_match(lab_names, c("creatinine", "^cr$"))

# Identify SpO2 vital name
vital_names <- unique(na.omit(vitals$vital_category))
spo2_name <- pick_first_match(vital_names, c("spo2", "o2 sat", "oxygen saturation"))

message("Detected names:",
        "\n  PaO2 lab_name: ", pao2_name,
        "\n  Creatinine lab_name: ", creat_name,
        "\n  SpO2 vital_name: ", spo2_name)

# Build FiO2 time series from respiratory_support (best available in your extract)
fio2_ts <- rs %>%
  transmute(hospitalization_id,
            t = safe_posix(recorded_dttm),
            fio2_frac = fio2_to_frac(fio2_set),
            device_category = as.character(device_category)) %>%
  filter(!is.na(t)) %>%
  arrange(hospitalization_id, t)

# PaO2 within 72h of ICU in, link nearest FiO2 (within 2h)
pao2_ts <- labs %>%
  filter(!is.na(pao2_name), lab_name == pao2_name) %>%
  transmute(hospitalization_id,
            t = coalesce(lab_collect_dttm, lab_result_dttm),
            pao2 = lab_value_num) %>%
  filter(!is.na(t), !is.na(pao2))


# SpO2 within 72h, link nearest FiO2 (within 2h)
spo2_ts <- vitals %>%
  filter(vital_category == "spo2") %>%
  transmute(hospitalization_id,
            t = recorded_dttm,
            spo2 = spo2_to_prop(vital_value)) %>%
  filter(!is.na(t), !is.na(spo2))

# function: nearest join within window per hospitalization
nearest_join <- function(x, y, by="hospitalization_id", max_diff_hours=2) {
  # x: events (has t), y: reference (has t)
  x <- x %>% arrange(.data[[by]], t)
  y <- y %>% arrange(.data[[by]], t)
  out <- x %>%
    left_join(y, by=by, suffix=c("", "_y")) %>%
    mutate(dt = abs(as.numeric(difftime(t_y, t, units="hours")))) %>%
    filter(!is.na(dt), dt <= max_diff_hours) %>%
    group_by(.data[[by]], t) %>%
    slice_min(dt, n=1, with_ties=FALSE) %>%
    ungroup() %>%
    select(-dt)
  out
}

# Restrict to 72h from ICU in (or you can do 7d post-op for ARF_7d)
t0_only <- t0_tbl %>% select(hospitalization_id, t0, proc_time)

pao2_72 <- pao2_ts %>%
  inner_join(t0_only, by="hospitalization_id") %>%
  mutate(dt_h = as.numeric(difftime(t, t0, units="hours"))) %>%
  filter(dt_h >= 0, dt_h <= 72) %>%
  select(hospitalization_id, t, pao2)

spo2_72 <- spo2_ts %>%
  inner_join(t0_only, by="hospitalization_id") %>%
  mutate(dt_h = as.numeric(difftime(t, t0, units="hours"))) %>%
  filter(dt_h >= 0, dt_h <= 72) %>%
  select(hospitalization_id, t, spo2)

fio2_72 <- fio2_ts %>%
  inner_join(t0_only, by="hospitalization_id") %>%
  mutate(dt_h = as.numeric(difftime(t, t0, units="hours"))) %>%
  filter(dt_h >= 0, dt_h <= 72) %>%
  select(hospitalization_id, t, fio2_frac, device_category)

pf <- nearest_join(pao2_72, fio2_72, max_diff_hours = 2) %>%
  mutate(pf_ratio = pao2 / fio2_frac)

sf <- nearest_join(spo2_72, fio2_72, max_diff_hours = 2) %>%
  mutate(sf_ratio = (spo2*100) / (fio2_frac*100))  # ratio in same units cancels; kept explicit

ox_summary <- full_join(
  pf %>% group_by(hospitalization_id) %>% summarise(worst_pf_72 = min(pf_ratio, na.rm=TRUE), .groups="drop"),
  sf %>% group_by(hospitalization_id) %>% summarise(worst_sf_72 = min(sf_ratio, na.rm=TRUE), .groups="drop"),
  by="hospitalization_id"
)

# ARDS-like physiology proxy within 72h:
# If PaO2 available: PF < 300; else SF < 315
ox_summary <- ox_summary %>%
  mutate(
    ards_proxy_72 = case_when(
      !is.na(worst_pf_72) ~ worst_pf_72 < 300,
      is.na(worst_pf_72) & !is.na(worst_sf_72) ~ worst_sf_72 < 315,
      TRUE ~ NA
    )
  )

# ---------------------------------------
# 4) ARF within 7 days post-op (support + FiO2 + hypoxemia)
# ---------------------------------------
rs_7d <- rs %>%
  inner_join(t0_tbl %>% select(hospitalization_id, proc_time), by="hospitalization_id") %>%
  mutate(dt_h = as.numeric(difftime(recorded_dttm, proc_time, units="hours"))) %>%
  filter(dt_h >= 0, dt_h <= 7*24) %>%
  mutate(
    fio2_frac = fio2_to_frac(fio2_set),
    state = map_state(device_category)
  ) %>%
  group_by(hospitalization_id) %>%
  summarise(
    imv_any_7d = any(state == "IMV", na.rm=TRUE),
    adv_any_7d = any(state %in% c("HFNC","NIV","IMV"), na.rm=TRUE),
    max_fio2_7d = max(fio2_frac, na.rm=TRUE),
    arf_support_7d = imv_any_7d | (adv_any_7d & max_fio2_7d >= 0.60),
    .groups="drop"
  ) %>%
  mutate(max_fio2_7d = ifelse(is.infinite(max_fio2_7d), NA_real_, max_fio2_7d))

# Hypoxemia component in 7d window (reuse oxygenation but expand to 7d if desired)
# For now: use worst 72h as early physiology; you can extend to 7d by mirroring code above.
phenos <- rs_7d %>%
  left_join(ox_summary, by="hospitalization_id") %>%
  mutate(
    arf_7d = case_when(
      arf_support_7d ~ TRUE,
      !is.na(ards_proxy_72) & ards_proxy_72 ~ TRUE,
      TRUE ~ FALSE
    )
  )

# ---------------------------------------
# 5) Ventilatory course outcomes (initiation, reintubation, escalation)
# ---------------------------------------
# IMV initiation after t0
imv_events <- rs %>%
  inner_join(t0_tbl %>% select(hospitalization_id, t0), by="hospitalization_id") %>%
  mutate(dt_h = as.numeric(difftime(recorded_dttm, t0, units="hours"))) %>%
  filter(dt_h >= 0, dt_h <= 7*24) %>%
  mutate(state = map_state(device_category)) %>%
  arrange(hospitalization_id, recorded_dttm) %>%
  group_by(hospitalization_id) %>%
  summarise(
    imv_initiated = any(state == "IMV", na.rm=TRUE),
    imv_start_time = ifelse(any(state=="IMV", na.rm=TRUE),
                            min(recorded_dttm[state=="IMV"], na.rm=TRUE),
                            as.POSIXct(NA)),
    .groups="drop"
  )

# Extubation failure / reintubation proxy:
# IMV -> non-IMV -> IMV within 48h of leaving IMV (requires sequences; do simple within 7d)
reintub <- rs %>%
  inner_join(t0_tbl %>% select(hospitalization_id, t0), by="hospitalization_id") %>%
  mutate(dt_h = as.numeric(difftime(recorded_dttm, t0, units="hours"))) %>%
  filter(dt_h >= 0, dt_h <= 7*24) %>%
  mutate(state = map_state(device_category)) %>%
  arrange(hospitalization_id, recorded_dttm) %>%
  group_by(hospitalization_id) %>%
  summarise(
    # crude: at least two distinct IMV episodes separated by >= 2h not on IMV
    reintubation = {
      imv_idx <- which(state == "IMV")
      if (length(imv_idx) < 2) FALSE else {
        gaps <- diff(recorded_dttm[imv_idx])
        any(as.numeric(gaps, units="hours") >= 2)
      }
    },
    .groups="drop"
  )

# Escalation within 72h: max severity in last 24h - first 24h
esc <- rs_seq %>%
  mutate(sev = as.integer(state)) %>%
  group_by(hospitalization_id) %>%
  summarise(
    max_0_24 = max(sev[bin <= 24], na.rm=TRUE),
    max_48_72 = max(sev[bin >= 48], na.rm=TRUE),
    escalation_72 = (max_48_72 - max_0_24) >= 2,
    .groups="drop"
  )

# ---------------------------------------
# 6) AKI within 7 days post-op (creatinine)
# ---------------------------------------
aki <- NULL
if (!is.na(creat_name)) {
  cr_ts <- labs %>%
    filter(lab_name == creat_name) %>%
    transmute(hospitalization_id,
              t = coalesce(lab_collect_dttm, lab_result_dttm),
              cr = lab_value_num) %>%
    filter(!is.na(t), !is.na(cr))
  
  # baseline: minimum in [-7d, 0] relative to proc_time; fallback: first post-op
  cr_bl <- cr_ts %>%
    inner_join(t0_tbl %>% select(hospitalization_id, proc_time), by="hospitalization_id") %>%
    mutate(dt_h = as.numeric(difftime(t, proc_time, units="hours"))) %>%
    group_by(hospitalization_id) %>%
    summarise(
      cr_baseline = {
        pre <- cr[dt_h >= -7*24 & dt_h <= 0]
        if (length(pre) > 0) min(pre, na.rm=TRUE) else {
          post <- cr[dt_h >= 0]
          if (length(post) > 0) post[which.min(dt_h[dt_h>=0])][1] else NA_real_
        }
      },
      .groups="drop"
    )
  
  cr_post <- cr_ts %>%
    inner_join(t0_tbl %>% select(hospitalization_id, proc_time), by="hospitalization_id") %>%
    mutate(dt_h = as.numeric(difftime(t, proc_time, units="hours"))) %>%
    filter(dt_h >= 0, dt_h <= 7*24)
  
  aki <- cr_post %>%
    left_join(cr_bl, by="hospitalization_id") %>%
    group_by(hospitalization_id) %>%
    summarise(
      cr_baseline = first(cr_baseline),
      cr_max_7d = max(cr, na.rm=TRUE),
      aki_7d = ifelse(!is.na(cr_baseline),
                      (cr_max_7d >= 1.5*cr_baseline),
                      NA),
      .groups="drop"
    )
}

# ---------------------------------------
# 7) Assemble analysis dataset
# ---------------------------------------

hosp_core <- clif_tables[["clif_hospitalization"]] %>%
  rename_with(tolower) %>%
  select(hospitalization_id, discharge_category)

dat <- cohort_lung_expo %>%
  left_join(rs_class, by="hospitalization_id") %>%
  left_join(phenos, by="hospitalization_id") %>%
  left_join(imv_events, by="hospitalization_id") %>%
  left_join(reintub, by="hospitalization_id") %>%
  left_join(esc, by="hospitalization_id")

dat <- dat %>%
  left_join(hosp_core, by = "hospitalization_id")

if (!is.null(aki)) dat <- dat %>% left_join(aki, by="hospitalization_id")

# define outcomes
dat <- dat %>%
  mutate(
    discharge_category = tolower(coalesce(discharge_category, "")),
    death_in_hosp = case_when(
      str_detect(discharge_category, "expired") ~ TRUE,
      str_detect(discharge_category, "death") ~ TRUE,
      str_detect(discharge_category, "hospice") ~ TRUE,
      TRUE ~ FALSE
    ),
    icu_los_days = icu_los_hours / 24
  )


# ---------------------------------------
# 8) Exposure association models
# ---------------------------------------
# Pick exposure variables (adjust to your names)
expo_vars <- c("pm25_1y","no2_1y","pm25_3y","no2_3y")
expo_vars <- expo_vars[expo_vars %in% names(dat)]

# Standardize exposures (per SD) for comparability
for (v in expo_vars) {
  dat[[paste0(v, "_z")]] <- as.numeric(scale(dat[[v]]))
}

# (A) Exposure -> Trajectory class (multinomial)
# Use class reference = "Stable low support"
if ("class" %in% names(dat)) {
  # minimal confounders already in cohort: age_years, sex_category, race_category
  # Add more later (comorbidity flags) if you derive them.
  f_mult <- as.formula("class ~ pm25_mean_lag + no2_mean_lag + age_years + sex_category + race_category")
  m_mult <- nnet::multinom(f_mult, data = dat, na.action = na.omit, trace = FALSE)
  
  mult_tidy <- broom::tidy(m_mult, exponentiate = TRUE, conf.int = TRUE)
  print(mult_tidy)
}

# (B) Exposure -> ARF_7d (logistic)
if ("arf_7d" %in% names(dat)) {
  m_arf <- glm(arf_7d ~ pm25_mean_lag + no2_mean_lag + age_years + sex_category,
               data = dat, family = binomial())
  
  print(broom::tidy(m_arf, exponentiate = TRUE, conf.int = TRUE))
}

# (C) Exposure -> ARDS proxy (logistic; among those with oxygenation)
if ("ards_proxy_72" %in% names(dat)) {
  m_ards <- glm(ards_proxy_72 ~ pm25_mean_lag + no2_mean_lag + age_years + sex_category,
                data = dat, family = binomial())
  print(broom::tidy(m_ards, exponentiate = TRUE, conf.int = TRUE))
}

# (D) Exposure -> ICU LOS (log-normal regression; robust/simple)
# (Alternative: Cox time-to-ICU-discharge if you build event times)
if ("icu_los_days" %in% names(dat)) {
  dat <- dat %>% mutate(log_icu_los = log1p(icu_los_days))
  m_los <- lm(log_icu_los ~ pm25_mean_lag + no2_mean_lag + age_years + sex_category,
              data = dat)
  print(broom::tidy(m_los, conf.int = TRUE))
}

# (E) Exposure -> In-hospital death (logistic)
if ("death_in_hosp" %in% names(dat)) {
  m_death <- glm(death_in_hosp ~ pm25_mean_lag + no2_mean_lag + age_years + sex_category,
                 data = dat, family = binomial())
  print(broom::tidy(m_death, exponentiate = TRUE, conf.int = TRUE))
}

# ---------------------------------------
# 9) Quick plots (exposure vs class)
# ---------------------------------------
if (all(c("class","pm25_mean_lag") %in% names(dat))) {
  ggplot(dat, aes(x = class, y = pm25_mean_lag)) +
    geom_boxplot() +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle=20, hjust=1)) +
    labs(title = "PM2.5 (1y pre-op) by respiratory trajectory class", x = NULL, y = "PM2.5")
}

if (all(c("class","no2_mean_lag") %in% names(dat))) {
  ggplot(dat, aes(x = class, y = no2_mean_lag)) +
    geom_boxplot() +
    theme_minimal(base_size = 14) +
    theme(axis.text.x = element_text(angle=20, hjust=1)) +
    labs(title = "PM2.5 (1y pre-op) by respiratory trajectory class", x = NULL, y = "PM2.5")
}


# ---------- Output directory ----------
sanitize_tag <- function(x) {
  x <- if (is.null(x)) "SITE" else as.character(x)
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "SITE" else x
}

SITE_NAME   <- sanitize_tag(site_name)
SYSTEM_DATE <- format(Sys.Date(), "%Y%m%d")
out_dir <- file.path(repo, "output", paste0("run_", SITE_NAME, "_", SYSTEM_DATE))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

save_plot <- function(p, filename, width = 10, height = 6, dpi = 320) {
  path <- file.path(out_dir, filename)
  ggsave(path, plot = p, width = width, height = height, dpi = dpi, units = "in")
  message("Saved: ", path)
}

save_csv <- function(df, filename) {
  path <- file.path(out_dir, filename)
  readr::write_csv(df, path)
  message("Saved: ", path)
}

pm_by_class <- dat %>%
  group_by(class) %>%
  summarise(
    n = sum(!is.na(pm25_mean_lag)),
    pm25_median = median(pm25_mean_lag, na.rm = TRUE),
    pm25_iqr = IQR(pm25_mean_lag, na.rm = TRUE),
    .groups = "drop"
  )

p_pm <- ggplot(dat, aes(x = class, y = pm25_mean_lag)) +
  geom_boxplot() +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle=20, hjust=1)) +
  labs(title = "PM2.5 (1y pre-op) by respiratory trajectory class", x = NULL, y = "PM2.5")

save_plot(p_pm, "fig_pm25_by_traj_class.png", width = 10, height = 5.5)
save_csv(pm_by_class, "table_pm25_by_traj_class_summary.csv")

no2_by_class <- dat %>%
  group_by(class) %>%
  summarise(
    n = sum(!is.na(no2_mean_lag)),
    pm25_median = median(no2_mean_lag, na.rm = TRUE),
    pm25_iqr = IQR(no2_mean_lag, na.rm = TRUE),
    .groups = "drop"
  )

p_no2 <- ggplot(dat, aes(x = class, y = no2_mean_lag)) +
  geom_boxplot() +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle=20, hjust=1)) +
  labs(title = "NO2 (1y pre-op) by respiratory trajectory class", x = NULL, y = "PM2.5")

save_plot(p_pm, "fig_no2_by_traj_class.png", width = 10, height = 5.5)
save_csv(pm_by_class, "table_no2_by_traj_class_summary.csv")


saveRDS(dat, file.path(out_dir, "analysis_dat_with_phenos_and_traj.rds"))
if (exists("m_arf")) saveRDS(m_arf, file.path(out_dir, "model_arf_logistic.rds"))
if (exists("m_mult")) saveRDS(m_mult, file.path(out_dir, "model_traj_multinom.rds"))
if (exists("m_los")) saveRDS(m_mult, file.path(out_dir, "model_traj_los.rds"))
if (exists("m_death")) saveRDS(m_mult, file.path(out_dir, "model_death.rds"))

m_arf

summary(m_arf)
summary(m_los)
summary(m_mult)
summary(m_death)

dat <- dat %>%
  left_join(rs_at_t0 %>% transmute(hospitalization_id, rs_t0 = rs_group), by="hospitalization_id") %>%
  mutate(imv_t0 = (rs_t0 == "IMV"))

rs2 <- rs %>%
  inner_join(t0_tbl %>% select(hospitalization_id, t0), by="hospitalization_id") %>%
  mutate(
    dt_h = as.numeric(difftime(recorded_dttm, t0, units="hours")),
    state = map_state(device_category),
    sev = match(state, c("None","RA","LowO2","HFNC","NIV","IMV","Other"))
  ) %>%
  filter(dt_h >= 0, dt_h <= 7*24) %>%
  arrange(hospitalization_id, recorded_dttm)

library(dplyr)
library(lubridate)

imv_course <- rs2 %>%
  arrange(hospitalization_id, recorded_dttm) %>%
  group_by(hospitalization_id) %>%
  group_modify(~{
    df <- .x
    
    # Scalars for this hospitalization
    imv_t0    <- any(df$state == "IMV" & df$dt_h <= 2, na.rm = TRUE)
    imv_any72 <- any(df$state == "IMV" & df$dt_h <= 72, na.rm = TRUE)
    
    # First IMV time (absolute + relative)
    if (any(df$state == "IMV" & df$dt_h <= 72, na.rm = TRUE)) {
      t_imv_abs <- df$recorded_dttm[df$state == "IMV" & df$dt_h <= 72][1]
      imv_start <- min(df$dt_h[df$state == "IMV" & df$dt_h <= 72], na.rm = TRUE)
    } else {
      t_imv_abs <- as.POSIXct(NA)
      imv_start <- NA_real_
    }
    
    # Extubation proxy: first non-IMV after first IMV
    if (!is.na(t_imv_abs)) {
      idx_ext <- which(df$recorded_dttm > t_imv_abs & df$state != "IMV")
      if (length(idx_ext) > 0) {
        t_ext_abs <- df$recorded_dttm[idx_ext[1]]
        extub_time <- as.numeric(difftime(t_ext_abs, df$t0[1], units = "hours"))
      } else {
        t_ext_abs <- as.POSIXct(NA)
        extub_time <- NA_real_
      }
    } else {
      t_ext_abs <- as.POSIXct(NA)
      extub_time <- NA_real_
    }
    
    extub_24h <- !is.na(extub_time) && extub_time <= 24
    extub_48h <- !is.na(extub_time) && extub_time <= 48
    
    # Reintubation within 48h of extubation
    reintub_48h <- FALSE
    if (!is.na(t_ext_abs)) {
      reintub_48h <- any(df$state == "IMV" &
                           df$recorded_dttm >  t_ext_abs &
                           df$recorded_dttm <= t_ext_abs + dhours(48),
                         na.rm = TRUE)
    }
    
    tibble(
      imv_t0 = imv_t0,
      imv_any_72 = imv_any72,
      imv_start_h = imv_start,
      extub_time_h = extub_time,
      extub_24h = extub_24h,
      extub_48h = extub_48h,
      reintub_48h = reintub_48h
    )
  }) %>%
  ungroup()


dat <- dat %>%
  left_join(imv_course %>%
              transmute(hospitalization_id,
                        imv_t0_flag = imv_t0,
                        imv_any_72 = imv_any_72,
                        extub_24h = extub_24h,
                        extub_48h = extub_48h,
                        reintub_48h = reintub_48h),
            by = "hospitalization_id")


dat <- dat %>%
  mutate(
    class_imv = case_when(
      imv_t0_flag & extub_24h ~ "IMV rapid extub (≤24h)",
      imv_t0_flag & !extub_24h & extub_48h ~ "IMV extub 24–48h",
      imv_t0_flag & reintub_48h ~ "IMV extub failure (reintub ≤48h)",
      imv_t0_flag ~ "IMV prolonged/unclear",
      TRUE ~ NA_character_
    ),
    class_nonimv = case_when(
      !imv_t0_flag & imv_any_72 ~ "Escalated to IMV",
      !imv_t0_flag & arf_7d ~ "ARF without IMV",
      !imv_t0_flag ~ "No IMV / stable",
      TRUE ~ NA_character_
    )
  )


dat_imv <- dat %>% filter(imv_t0)

m_extub24 <- glm(extub_24h ~ pm25_mean_lag + no2_mean_lag + age_years + sex_category,
                 data = dat_imv, family = binomial())
saveRDS(m_extub24, file.path(out_dir, "model_extub24_imv_t0.rds"))

dat_non <- dat %>% filter(!imv_t0)

m_escal <- glm(imv_any_72 ~ pm25_mean_lag + no2_mean_lag + age_years + sex_category,
               data = dat_non, family = binomial())
saveRDS(m_escal, file.path(out_dir, "model_escalation_to_imv_nonimv_t0.rds"))

summary(m_extub24)
summary(m_escal)

dat$pm25_z <- scale(dat$pm25_mean_lag)
dat$no2_z  <- scale(dat$no2_mean_lag)

dat <- dat %>% left_join(cohort_lung %>% select(hospitalization_id, surgery_category), by="hospitalization_id")

# stratified
dat_non <- dat %>% filter(!imv_t0_flag)
by(dat_non, dat_non$surgery_category, function(df) {
  glm(imv_any_72 ~ pm25_mean_lag + no2_mean_lag + age_years + sex_category,
      data=df, family=binomial())
})















