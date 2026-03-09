# ================================================================================================
# ICU Post–Lung Cancer Resection Cohort Builder (CLIF) | PI: Peter Graffy
#
# Rebuilt cohort logic (robust to site variation in procedure_code_format):
#   1) Identify lung cancer dx POA=1 from clif_hospital_diagnosis (ICD9/10 prefix match)
#   2) Identify lung resection procedures from clif_patient_procedures by matching:
#        - CPT numeric codes to curated CPT list
#        - ICD10PCS alphanumeric codes to curated PCS list
#      (Does NOT rely on procedure_code_format to infer code system.)
#   3) Determine who went to ICU (clif_adt location_category contains "ICU")
#   4) Link procedure -> ICU: keep closest ICU admission AFTER procedure within window
#
# Outputs (repo/output/run_[SITE]_[DATE]/):
#   cohort_lung_resection_icu.csv
#   exclusion_lung_resection_icu.csv
#   flow_lung_resection_icu.csv
# ================================================================================================

suppressPackageStartupMessages({
  library(fst)
  library(here)
  library(tidyverse)
  library(arrow)
  library(stringr)
  library(lubridate)
  library(readr)
  library(glue)
  library(arrow)
  library(scales)
})

# ---------- Project config ----------
source("utils/config.R")
stopifnot(exists("config"))

repo        <- config$repo
site_name   <- config$site_name
tables_path <- config$tables_path
file_type   <- config$file_type

stopifnot(!is.null(repo), nzchar(repo))
stopifnot(!is.null(tables_path), nzchar(tables_path))

cat("Site Name:", site_name, "\n")
cat("Tables Path:", tables_path, "\n")
cat("File Type:", file_type, "\n")

tables_path <- normalizePath(tables_path, mustWork = TRUE)

# ---------- Parameters ----------
START_DATE <- as.POSIXct("2018-01-01 00:00:00", tz = "UTC")
END_DATE   <- as.POSIXct("2024-12-31 23:59:59", tz = "UTC")

ADULT_AGE_YEARS <- 18
PROC_TO_ICU_MAX_H <- 48
T0_ANCHOR <- "icu_in" # "icu_in" or "proc"

CODES_PATH <- file.path(repo, "outlier-thresholds", "table.tsv")
stopifnot(file.exists(CODES_PATH))

# ---------- Helpers ----------
safe_posix <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  if (is.numeric(x)) return(as.POSIXct(x, origin = "1970-01-01", tz = "UTC"))
  suppressWarnings(as.POSIXct(x, tz = "UTC"))
}
safe_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  suppressWarnings(as.Date(x))
}
sanitize_tag <- function(x) {
  x <- if (is.null(x)) "SITE" else as.character(x)
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "SITE" else x
}
norm_code <- function(x) {
  x <- toupper(as.character(x))
  str_replace_all(x, "[^A-Z0-9]", "")
}
code_matches_any_prefix <- function(code_vec, prefixes) {
  code_vec <- norm_code(code_vec)
  prefixes <- norm_code(prefixes)
  vapply(code_vec, function(cd) any(str_starts(cd, prefixes)), logical(1))
}

# ---------- Output folder ----------
SITE_NAME   <- sanitize_tag(site_name)
SYSTEM_DATE <- format(Sys.Date(), "%Y%m%d")

out_dir <- file.path(repo, "output", paste0("run_", SITE_NAME, "_", SYSTEM_DATE))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

make_name <- function(stub, ext="csv") paste0(stub, "_", SITE_NAME, "_", SYSTEM_DATE, ".", ext)

# ---------- File discovery / loading ----------
exts <- strsplit(file_type, "[/|,; ]+")[[1]]
exts <- exts[nzchar(exts)]
if (length(exts) == 0) exts <- c("csv","parquet","fst")
ext_pat <- paste0("\\.(", paste(unique(exts), collapse = "|"), ")$")

all_files <- list.files(
  path = tables_path,
  pattern = ext_pat,
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(all_files) == 0) {
  stop("No files with extensions {", paste(exts, collapse = ", "), "} were found under: ", tables_path)
}

bn <- basename(all_files)
looks_clif <- grepl("^clif_.*", bn, ignore.case = TRUE)
base_no_ext <- tools::file_path_sans_ext(tolower(bn))
base_norm <- ifelse(looks_clif, base_no_ext, paste0("clif_", base_no_ext))
found_map <- stats::setNames(all_files, base_norm)

required_raw <- c("patient","hospitalization","adt","hospital_diagnosis","patient_procedures", "respiratory_support", "labs", "vitals")
required_files <- paste0("clif_", required_raw)
missing <- setdiff(required_files, names(found_map))
if (length(missing) > 0) {
  cat("Detected CLIF-like files:\n"); print(sort(unique(names(found_map))))
  stop("Missing required tables: ", paste(missing, collapse = ", "))
}

clif_paths <- found_map[required_files]

read_any <- function(path) {
  ext <- tolower(tools::file_ext(path))
  switch(ext,
         "csv"     = readr::read_csv(path, show_col_types = FALSE),
         "parquet" = arrow::read_parquet(path),
         "fst"     = fst::read_fst(path, as.data.table = FALSE),
         stop("Unsupported extension: ", ext))
}

clif_tables <- lapply(clif_paths, read_any)
names(clif_tables) <- required_files
cat("Loaded tables: ", paste(names(clif_tables), collapse = ", "), "\n")

get_min <- function(tbl_name, cols) {
  nm <- paste0("clif_", tbl_name)
  out <- clif_tables[[nm]] %>% rename_with(tolower)
  cols_keep <- intersect(tolower(cols), names(out))
  out %>% select(any_of(cols_keep))
}

# ---------- Load codes table (dx prefixes) ----------
codes <- readr::read_tsv(CODES_PATH, show_col_types = FALSE) %>%
  rename_with(tolower)

stopifnot(all(c("code type","code") %in% names(codes)))

codes <- codes %>%
  mutate(
    code_system = toupper(str_replace_all(`code type`, "[^A-Z0-9]", "")),
    code = norm_code(code)
  )

lung_dx_prefixes <- codes %>%
  filter(code_system %in% c("ICD10CM","ICD9CM")) %>%
  pull(code) %>%
  unique()

cat("Loaded lung dx prefix count: ", length(lung_dx_prefixes), "\n", sep="")

# ---------- Curated resection procedure codes (DO NOT infer from code_system columns) ----------
resection_cpt <- norm_code(c(
  "32480","32482","32663","32484","32505","32666","32669","32440","32442"
))
resection_icd10pcs <- norm_code(c(
  "0BTC0ZZ","0BTC4ZZ","0BTD0ZZ","0BTD4ZZ","0BTG0ZZ","0BTG4ZZ",
  "0BTH0ZZ","0BTH4ZZ","0BTJ0ZZ","0BTJ4ZZ","0BTK0ZZ","0BTK4ZZ"
))

# surgery categories (from code patterns)
surg_from_code <- function(proc_code) {
  pc <- norm_code(proc_code)
  case_when(
    pc %in% c("32440","32442") ~ "Pneumonectomy",
    pc %in% c("32480","32482","32663") ~ "Lobectomy",
    pc %in% c("32484","32669") ~ "Segmentectomy",
    pc %in% c("32505","32666") ~ "Wedge",
    str_starts(pc, "0BTJ") ~ "Pneumonectomy",
    str_starts(pc, "0BTK") ~ "Pneumonectomy",
    str_starts(pc, "0BTC") ~ "Lobectomy",
    str_starts(pc, "0BTD") ~ "Lobectomy",
    str_starts(pc, "0BTG") ~ "Lobectomy",
    str_starts(pc, "0BTH") ~ "Lobectomy",
    TRUE ~ "Other"
  )
}

# ---------- Minimal tables ----------
patient <- get_min("patient",
                   c("patient_id","birth_date","sex_category","race_category","ethnicity_category")) %>%
  mutate(birth_date = safe_date(birth_date))

hospitalization <- get_min("hospitalization",
                           c("patient_id",
                             "hospitalization_id",
                             "admission_dttm",
                             "discharge_dttm",
                             "discharge_category",
                             "age_at_admission",
                             "zipcode_nine_digit",
                             "zipcode_five_digit",
                             "census_tract",
                             "county_code")) %>%
  mutate(
    admission_dttm = safe_posix(admission_dttm),
    discharge_dttm = safe_posix(discharge_dttm),
    discharge_category = as.character(discharge_category)
  )

hospitalization <- hospitalization %>%
  mutate(
    discharge_l = tolower(coalesce(discharge_category, "")),
    
    death_in_hosp = str_detect(discharge_l, "expired|death"),
    
    hospice_discharge = str_detect(discharge_l, "hospice"),
    
    death_or_hospice = death_in_hosp | hospice_discharge
  )


adt <- get_min("adt",
               c("hospitalization_id","in_dttm","out_dttm","location_category","location_type")) %>%
  mutate(
    in_dttm  = safe_posix(in_dttm),
    out_dttm = safe_posix(out_dttm),
    location_category = as.character(location_category)
  )

hospital_dx <- get_min("hospital_diagnosis",
                       c("hospitalization_id","diagnosis_code","diagnosis_code_format","poa_present","diagnosis_primary")) %>%
  mutate(
    icd_code = norm_code(diagnosis_code),
    poa_present = suppressWarnings(as.integer(poa_present))
  )

patient_proc <- get_min("patient_procedures",
                        c("hospitalization_id","procedure_code","procedure_code_format","procedure_billed_dttm")) %>%
  mutate(
    procedure_code = norm_code(procedure_code),
    procedure_billed_dttm = safe_posix(procedure_billed_dttm),
    proc_time = procedure_billed_dttm
  )

# ---------- Step A: base hospitalizations in time window + adults + geo ----------
base <- hospitalization %>%
  filter(!is.na(admission_dttm),
         admission_dttm >= START_DATE,
         admission_dttm <= END_DATE) %>%
  left_join(patient %>% select(patient_id, birth_date, sex_category, race_category, ethnicity_category),
            by = "patient_id") %>%
  mutate(
    age_years = coalesce(
      suppressWarnings(as.numeric(age_at_admission)),
      ifelse(!is.na(birth_date),
             as.numeric(floor((as.Date(admission_dttm) - birth_date)/365.25)), NA_real_)
    ),
    has_demo = !(is.na(age_years) | is.na(sex_category) | is.na(race_category)),
    adult    = !is.na(age_years) & age_years >= ADULT_AGE_YEARS,
    has_geo  = !is.na(census_tract) | !is.na(zipcode_nine_digit) | !is.na(zipcode_five_digit) | !is.na(county_code)
  )

# ---------- Step B: lung cancer dx POA=1 ----------
lung_poa <- hospital_dx %>%
  filter(poa_present == 1) %>%
  filter(code_matches_any_prefix(icd_code, lung_dx_prefixes)) %>%
  distinct(hospitalization_id) %>%
  mutate(has_lung_cancer_poa = TRUE)

# ---------- Step C: resection procedures (match by code shape, not code_system) ----------
# CPT codes are numeric; PCS codes are alphanumeric and here start with 0B...
proc_resection <- patient_proc %>%
  filter(!is.na(proc_time), !is.na(procedure_code), procedure_code != "") %>%
  mutate(
    is_cpt = str_detect(procedure_code, "^[0-9]+$"),
    is_pcs = str_detect(procedure_code, "^[A-Z0-9]+$") & str_starts(procedure_code, "0B")
  ) %>%
  filter(
    (is_cpt & procedure_code %in% resection_cpt) |
      (is_pcs & procedure_code %in% resection_icd10pcs)
  ) %>%
  mutate(
    surgery_category = surg_from_code(procedure_code)
  ) %>%
  transmute(
    hospitalization_id,
    proc_time,
    procedure_code,
    surgery_category
  )

# ---------- Step D: ICU segments + bounds ----------
icu_segments <- adt %>%
  mutate(is_icu = str_detect(tolower(coalesce(location_category, "")), "icu")) %>%
  filter(is_icu, !is.na(in_dttm)) %>%
  transmute(
    hospitalization_id,
    icu_in_time = in_dttm,
    icu_out_time = out_dttm
  )

icu_bounds <- icu_segments %>%
  group_by(hospitalization_id) %>%
  summarize(
    first_icu_in = suppressWarnings(min(icu_in_time, na.rm = TRUE)),
    last_icu_out = suppressWarnings(max(icu_out_time, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    first_icu_in = ifelse(is.infinite(first_icu_in), NA, first_icu_in),
    last_icu_out = ifelse(is.infinite(last_icu_out), NA, last_icu_out),
    icu_los_hours = as.numeric(difftime(last_icu_out, first_icu_in, units = "hours"))
  )

# ---------- Step E: restrict to hospitalizations with dx + procedure ----------
eligible <- base %>%
  inner_join(lung_poa, by = "hospitalization_id") %>%
  inner_join(proc_resection, by = "hospitalization_id")

# ---------- Step F: link procedure to ICU (closest ICU admission after procedure) ----------
proc_icu_pairs <- eligible %>%
  inner_join(icu_segments, by = "hospitalization_id") %>%
  mutate(proc_to_icu_h = as.numeric(difftime(icu_in_time, proc_time, units = "hours"))) %>%
  filter(proc_to_icu_h >= 0) %>%
  filter(is.infinite(PROC_TO_ICU_MAX_H) | proc_to_icu_h <= PROC_TO_ICU_MAX_H) %>%
  group_by(hospitalization_id) %>%
  slice_min(proc_to_icu_h, n = 1, with_ties = FALSE) %>%
  ungroup()

# ---------- Final cohort ----------
cohort_lung <- proc_icu_pairs %>%
  left_join(icu_bounds, by = "hospitalization_id") %>%
  mutate(
    t0 = case_when(
      T0_ANCHOR == "icu_in" ~ icu_in_time,
      T0_ANCHOR == "proc"   ~ proc_time,
      TRUE ~ icu_in_time
    )
  ) %>%
  transmute(
    cohort = "lung_resection_icu",
    patient_id, hospitalization_id, adult, has_demo, has_geo,
    admission_dttm, discharge_dttm,
    first_icu_in, last_icu_out, icu_los_hours,
    age_years, sex_category, race_category, ethnicity_category,
    census_tract, county_code, zipcode_five_digit, zipcode_nine_digit,
    proc_time, procedure_code, surgery_category,
    icu_in_time, icu_out_time, proc_to_icu_h,
    t0
  ) %>%
  # apply final inclusion gates
  mutate(
    include = coalesce(adult, FALSE) & coalesce(has_demo, FALSE) & coalesce(has_geo, FALSE)
  ) %>%
  filter(include) %>%
  select(-include)

# ---------- Exclusions + flow ----------
# Start from base and track steps explicitly
step0 <- base
step1 <- step0 %>% filter(adult)
step2 <- step1 %>% filter(has_demo)
step3 <- step2 %>% filter(has_geo)
step4 <- step3 %>% inner_join(lung_poa, by="hospitalization_id")
step5 <- step4 %>% inner_join(proc_resection %>% distinct(hospitalization_id), by="hospitalization_id")
step6 <- step5 %>%
  inner_join(icu_segments %>% distinct(hospitalization_id), by="hospitalization_id")

flow_lung <- tibble(
  step = c(
    "Hospitalizations in date window",
    glue(">= {ADULT_AGE_YEARS} years"),
    "Demographics present",
    "Geography present",
    "Lung cancer POA present",
    "Resection procedure present",
    "Any ICU segment present"
  ),
  remaining = c(nrow(step0), nrow(step1), nrow(step2), nrow(step3), nrow(step4), nrow(step5), nrow(step6))
) %>%
  mutate(excluded_at_step = lag(remaining, default = remaining[1]) - remaining)

# Exclusions (first failing reason)
exclusion_lung <- base %>%
  mutate(
    reason = case_when(
      !adult ~ "Under 18 or missing age",
      !has_demo ~ "Missing demographics",
      !has_geo ~ "Missing geo code",
      TRUE ~ NA_character_
    )
  ) %>%
  # for those who pass demo/geo, check dx/proc/icu failure
  left_join(lung_poa %>% transmute(hospitalization_id, dx_ok = TRUE), by="hospitalization_id") %>%
  left_join(proc_resection %>% distinct(hospitalization_id) %>% mutate(proc_ok = TRUE), by="hospitalization_id") %>%
  left_join(icu_segments %>% distinct(hospitalization_id) %>% mutate(icu_ok = TRUE), by="hospitalization_id") %>%
  mutate(
    reason = coalesce(
      reason,
      ifelse(is.na(dx_ok), "No lung cancer dx present-on-admission (poa_present != 1 or no matching dx code)", NA_character_),
      ifelse(is.na(proc_ok), "No qualifying lung resection procedure code", NA_character_),
      ifelse(is.na(icu_ok), "No ICU stay (no ADT segment with ICU location_category)", NA_character_),
      "Other"
    )
  ) %>%
  filter(!(hospitalization_id %in% cohort_lung$hospitalization_id)) %>%
  select(patient_id, hospitalization_id, reason)

# ---------- Save outputs ----------
readr::write_csv(cohort_lung, file.path(out_dir, make_name("cohort_lung_resection_icu")))
readr::write_csv(exclusion_lung, file.path(out_dir, make_name("exclusion_lung_resection_icu")))
readr::write_csv(flow_lung, file.path(out_dir, make_name("flow_lung_resection_icu")))

cat("\nSaved outputs to: ", out_dir, "\n", sep="")
cat("  - ", make_name("cohort_lung_resection_icu"), "\n", sep="")
cat("  - ", make_name("exclusion_lung_resection_icu"), "\n", sep="")
cat("  - ", make_name("flow_lung_resection_icu"), "\n", sep="")

cat("\nCohort selection summary:\n")
cat("  Base hospitalizations:         ", nrow(step0), "\n", sep="")
cat("  Eligible (dx+proc):            ", nrow(step5), "\n", sep="")
cat("  Eligible + ICU:                ", nrow(step6), "\n", sep="")
cat("  Final cohort included:         ", nrow(cohort_lung), "\n", sep="")

# ---- configure paths ----
repo <- config$repo
SITE_NAME <- gsub("[^A-Za-z0-9]+", "_", config$site_name %||% "SITE")
SYSTEM_DATE <- format(Sys.Date(), "%Y%m%d")
out_dir <- file.path(repo, "output", paste0("analysis_", SITE_NAME, "_", SYSTEM_DATE))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

save_csv <- function(df, nm) readr::write_csv(df, file.path(out_dir, paste0(nm, ".csv")))
save_rds <- function(obj, nm) saveRDS(obj, file.path(out_dir, paste0(nm, ".rds")))

theme_pub <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          axis.title = element_text(face="bold"),
          plot.title = element_text(face="bold"))
}
save_plot <- function(p, filename, w=9, h=5.5, dpi=320) {
  ggsave(file.path(out_dir, filename), plot=p, width=w, height=h, dpi=dpi, units="in")
}

# ---- utility: robust county normalization ----
norm_county <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "[^0-9]", "")
  x <- ifelse(nchar(x) == 4, paste0("0", x), x)
  x
}

# ---- utility: safe numeric ----
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

expo_dir <- file.path(repo, "exposome")

# You may have csv/parquet/fst; adjust as needed
pm25 <- readr::read_csv(file.path(expo_dir, "pm25_county_year.csv"), show_col_types = FALSE)
no2  <- readr::read_csv(file.path(expo_dir, "no2_county_year.csv"),  show_col_types = FALSE)

# ---- REQUIRED: rename to county + year + value ----
# Adjust these three lines based on your actual columns
pm25 <- pm25 %>%
  rename(county = GEOID, year = year, pm25 = pm25_mean) %>%
  mutate(county = norm_county(county), year = as.integer(year), pm25 = as_num(pm25))

no2 <- no2 %>%
  rename(county = GEOID, year = year, no2 = no2_mean) %>%
  mutate(county = norm_county(county), year = as.integer(year), no2 = as_num(no2))

cohort <- cohort_lung %>%
  mutate(
    county = norm_county(county_code),
    proc_date = as.Date(proc_time),
    proc_year = as.integer(format(proc_date, "%Y"))
  )

stopifnot("county" %in% names(cohort), "proc_year" %in% names(cohort))

# helper: compute rolling mean from annual table
get_window_mean <- function(expo_df, value_col, years_back) {
  value_col <- rlang::ensym(value_col)
  expo_df %>%
    group_by(county) %>%
    arrange(year) %>%
    ungroup() %>%
    select(county, year, !!value_col)
}

pm25a <- get_window_mean(pm25, pm25, years_back = 3)
no2a  <- get_window_mean(no2,  no2,  years_back = 3)

# Build "target years" per hospitalization for 1y and 3y
cohort_expo <- cohort %>%
  mutate(
    y1 = proc_year - 1,
    y3_start = proc_year - 3,
    y3_end   = proc_year - 1
  )

# 1-year: match year == proc_year-1
pm25_1y <- cohort_expo %>%
  select(hospitalization_id, county, y1) %>%
  left_join(pm25a, by = c("county"="county", "y1"="year")) %>%
  rename(pm25_1y = pm25) %>%
  select(hospitalization_id, pm25_1y)

no2_1y <- cohort_expo %>%
  select(hospitalization_id, county, y1) %>%
  left_join(no2a, by = c("county"="county", "y1"="year")) %>%
  rename(no2_1y = no2) %>%
  select(hospitalization_id, no2_1y)

# 3-year: average years [proc_year-3, proc_year-1]
pm25_3y <- cohort_expo %>%
  select(hospitalization_id, county, y3_start, y3_end) %>%
  left_join(pm25a, by="county") %>%
  filter(year >= y3_start, year <= y3_end) %>%
  group_by(hospitalization_id) %>%
  summarise(pm25_3y = mean(pm25, na.rm=TRUE), .groups="drop")

no2_3y <- cohort_expo %>%
  select(hospitalization_id, county, y3_start, y3_end) %>%
  left_join(no2a, by="county") %>%
  filter(year >= y3_start, year <= y3_end) %>%
  group_by(hospitalization_id) %>%
  summarise(no2_3y = mean(no2, na.rm=TRUE), .groups="drop")

dat <- cohort %>%
  left_join(pm25_1y, by="hospitalization_id") %>%
  left_join(no2_1y,  by="hospitalization_id") %>%
  left_join(pm25_3y, by="hospitalization_id") %>%
  left_join(no2_3y,  by="hospitalization_id")

qc_expo <- dat %>%
  summarise(
    n = n(),
    pm25_1y_missing = mean(is.na(pm25_1y)),
    no2_1y_missing  = mean(is.na(no2_1y)),
    pm25_3y_missing = mean(is.na(pm25_3y)),
    no2_3y_missing  = mean(is.na(no2_3y))
  )

print(qc_expo)
save_csv(qc_expo, "qc_exposure_linkage")
save_csv(dat, "analysis_dat_with_exposures")

rs_raw <- clif_tables[["clif_respiratory_support"]] %>%
  rename_with(tolower) %>%
  mutate(
    recorded_dttm = safe_posix(recorded_dttm),
    device_category = as.character(device_category)
  )

map_rs_group <- function(device_category) {
  d <- as.character(device_category)
  case_when(
    d == "IMV" ~ "IMV",
    d %in% c("NIPPV","CPAP") ~ "NIV/CPAP",
    d == "High Flow NC" ~ "HFNC",
    d %in% c("NC","Nasal Cannula","Simple Mask","Non-Rebreather","Face Mask","Venturi") ~ "Low-flow O2",
    d %in% c("Room Air","RA") ~ "Room Air",
    is.na(d) | d == "" ~ "No RS recorded",
    TRUE ~ "Other"
  )
}

t0_tbl <- dat %>% transmute(hospitalization_id, t0 = icu_in_time)

rs2 <- rs_raw %>%
  inner_join(t0_tbl, by="hospitalization_id") %>%
  filter(!is.na(recorded_dttm), !is.na(t0)) %>%
  mutate(
    dt_h = as.numeric(difftime(recorded_dttm, t0, units="hours")),
    state = map_rs_group(device_category)
  ) %>%
  filter(dt_h >= -2, dt_h <= 72) %>%  # include a small pre-window for t0 assignment
  arrange(hospitalization_id, recorded_dttm)

save_rds(rs2, "rs_panel_72h")

severity <- function(state) {
  case_when(
    state == "IMV" ~ 5L,
    state == "NIV/CPAP" ~ 4L,
    state == "HFNC" ~ 3L,
    state == "Low-flow O2" ~ 2L,
    state == "Room Air" ~ 1L,
    TRUE ~ 0L
  )
}

imv_course <- rs2 %>%
  group_by(hospitalization_id) %>%
  summarise(
    imv_t0 = any(state == "IMV" & abs(dt_h) <= 2, na.rm=TRUE),
    imv_any_72 = any(state == "IMV" & dt_h >= 0 & dt_h <= 72, na.rm=TRUE),
    
    # time of first IMV (hours from t0), within 0–72h
    imv_start_h = if (any(state=="IMV" & dt_h>=0 & dt_h<=72, na.rm=TRUE)) {
      min(dt_h[state=="IMV" & dt_h>=0 & dt_h<=72], na.rm=TRUE)
    } else NA_real_,
    
    # extubation time: first non-IMV recorded after first IMV
    extub_h = {
      if (!any(state=="IMV" & dt_h>=0, na.rm=TRUE)) NA_real_ else {
        idx_imv <- which(state=="IMV" & dt_h>=0)
        first_imv_idx <- idx_imv[1]
        idx_post <- which(seq_along(state) > first_imv_idx & dt_h >= 0 & state != "IMV")
        if (length(idx_post)==0) NA_real_ else dt_h[idx_post[1]]
      }
    },
    
    extub_24h = !is.na(extub_h) & extub_h <= 24,
    extub_48h = !is.na(extub_h) & extub_h <= 48,
    
    # reintubation within 48h after extubation
    reintub_48h = {
      if (is.na(extub_h)) FALSE else {
        any(state=="IMV" & dt_h > extub_h & dt_h <= extub_h + 48, na.rm=TRUE)
      }
    },
    
    # max severity within 24h
    max_sev_24h = {
      d24 <- dt_h >= 0 & dt_h <= 24
      if (!any(d24, na.rm=TRUE)) NA_integer_ else max(severity(state[d24]), na.rm=TRUE)
    },
    
    .groups="drop"
  ) %>%
  mutate(
    max_support_24h = case_when(
      max_sev_24h == 5L ~ "IMV",
      max_sev_24h == 4L ~ "NIV/CPAP",
      max_sev_24h == 3L ~ "HFNC",
      max_sev_24h == 2L ~ "Low-flow O2",
      max_sev_24h == 1L ~ "Room Air",
      max_sev_24h == 0L ~ "None/Other",
      TRUE ~ NA_character_
    ),
    max_support_24h = factor(max_support_24h,
                             levels=c("Room Air","Low-flow O2","HFNC","NIV/CPAP","IMV","None/Other"))
  )

dat2 <- dat %>%
  left_join(imv_course, by="hospitalization_id")

rs_long <- rs_raw %>%
  inner_join(
    cohort_lung %>%
      select(hospitalization_id, icu_in_time, icu_out_time),
    by = "hospitalization_id"
  ) %>%
  mutate(
    recorded_dttm = safe_posix(recorded_dttm)
  ) %>%
  filter(
    recorded_dttm >= icu_in_time,
    recorded_dttm <= icu_out_time
  ) %>%
  arrange(hospitalization_id, recorded_dttm) %>%
  group_by(hospitalization_id) %>%
  mutate(
    next_time = lead(recorded_dttm),
    interval_hours = as.numeric(difftime(next_time, recorded_dttm, units="hours")),
    interval_hours = ifelse(is.na(interval_hours), 0, interval_hours)
  ) %>%
  ungroup()

imv_hours_total <- rs_long %>%
  filter(device_category == "IMV") %>%
  group_by(hospitalization_id) %>%
  summarise(
    imv_hours_total = sum(interval_hours, na.rm=TRUE),
    .groups="drop"
  )

rs_long_72 <- rs_long %>%
  mutate(
    hours_from_icu = as.numeric(difftime(recorded_dttm, icu_in_time, units="hours"))
  ) %>%
  filter(hours_from_icu <= 72)

imv_hours_72 <- rs_long_72 %>%
  filter(device_category == "IMV") %>%
  group_by(hospitalization_id) %>%
  summarise(
    imv_hours_72 = sum(interval_hours, na.rm=TRUE),
    .groups="drop"
  )

save_csv(imv_course, "resp_course_features")
save_csv(dat2, "analysis_dat_with_course_features")

# exposure scaling: IQR
iqr_scale <- function(x) {
  i <- IQR(x, na.rm=TRUE)
  (x - median(x, na.rm=TRUE)) / i
}

dat2 <- dat2 %>%
  mutate(
    pm25_iqr1y = iqr_scale(pm25_1y),
    no2_iqr1y  = iqr_scale(no2_1y),
    pm25_iqr3y = iqr_scale(pm25_3y),
    no2_iqr3y  = iqr_scale(no2_3y)
  )

dat_imv <- dat2 %>% filter(imv_t0)

dat_non <- dat2 %>% filter(!imv_t0)

m_esc_pm25 <- glm(imv_any_72 ~ pm25_iqr1y + age_years + sex_category + surgery_category,
                  data = dat_non, family = binomial())

m_esc_no2 <- glm(imv_any_72 ~ no2_iqr1y + age_years + sex_category + surgery_category,
                 data = dat_non, family = binomial())

save_rds(m_esc_pm25, "model_escalation_pm25_glm")
save_rds(m_esc_no2,  "model_escalation_no2_glm")
summary(m_esc_pm25)
summary(m_esc_no2)

p_expo <- dat2 %>%
  pivot_longer(c(pm25_1y, no2_1y), names_to="expo", values_to="value") %>%
  ggplot(aes(x=expo, y=value)) +
  geom_boxplot() +
  labs(title="Pre-op county-level exposures (1y)", x=NULL, y=NULL) +
  theme_pub()

save_plot(p_expo, "fig_exposure_boxplot.png", w=7, h=5)


dat2 <- dat2 %>%
  mutate(
    course = case_when(
      imv_t0 & extub_24h ~ "IMV rapid extub ≤24h",
      imv_t0 & !extub_24h & extub_48h ~ "IMV extub 24–48h",
      imv_t0 & !extub_48h ~ "IMV prolonged/unclear",
      !imv_t0 & imv_any_72 ~ "Escalated to IMV (≤72h)",
      !imv_t0 & !imv_any_72 ~ "No IMV in 72h",
      TRUE ~ NA_character_
    )
  )

p_pm25_course <- ggplot(dat2, aes(x=course, y=pm25_1y)) +
  geom_boxplot() +
  labs(title="PM2.5 (1y pre-op) by respiratory course", x=NULL, y="PM2.5") +
  theme_pub() +
  theme(axis.text.x = element_text(angle=25, hjust=1))

save_plot(p_pm25_course, "fig_pm25_by_course.png", w=10, h=5.5)

dat2 <- dat2 %>%
  left_join(
    hospitalization %>%
      select(hospitalization_id,
             discharge_category,
             death_in_hosp,
             hospice_discharge,
             death_or_hospice),
    by = "hospitalization_id"
  )


dat2 %>%
  summarise(
    n = n(),
    imv_t0 = sum(imv_t0, na.rm=TRUE),
    non_imv_t0 = sum(!imv_t0, na.rm=TRUE),
    extub_48h = sum(extub_48h, na.rm=TRUE),
    imv_any_72 = sum(imv_any_72, na.rm=TRUE),
    reintub_48h = sum(reintub_48h, na.rm=TRUE),
    deaths = sum(death_in_hosp, na.rm=TRUE)
  )

dat2 <- dat2 %>%
  mutate(
    hosp_los_days = as.numeric(difftime(discharge_dttm, admission_dttm, units="days"))
  )

dat2 <- dat2 %>%
  mutate(
    pm25_tertile = ntile(pm25_1y, 3),
    no2_tertile  = ntile(no2_1y, 3)
  )

tab_extub <- dat2 %>%
  filter(imv_t0) %>%
  group_by(pm25_tertile) %>%
  summarise(
    n = n(),
    extub_48h_rate = mean(extub_48h, na.rm=TRUE)
  )
tab_extub

p1 <- ggplot(tab_extub, aes(x=factor(pm25_tertile), y=extub_48h_rate)) +
  geom_col() +
  scale_y_continuous(labels=scales::percent_format()) +
  labs(title="Extubation ≤48h by PM2.5 tertile",
       x="PM2.5 tertile",
       y="Proportion extubated ≤48h") +
  theme_minimal()

p1



tab_matrix <- dat2 %>%
  filter(imv_t0) %>%
  select(pm25_tertile, extub_48h) %>%
  table()

fisher.test(tab_matrix)

dat2 %>%
  filter(!is.na(pm25_1y)) %>%
  ggplot(aes(x=factor(pm25_tertile), y=pm25_1y)) +
  geom_boxplot() +
  geom_jitter(width=.1, alpha=.6) +
  theme_minimal()

dat2 %>%
  mutate(pm25_tertile = ntile(pm25_1y, 3)) %>%
  count(pm25_tertile, death_in_hosp)

fisher.test(table(dat2$pm25_tertile, dat2$death_in_hosp))
fisher.test(table(dat2$no2_tertile, dat2$death_in_hosp))

table(dat2$pm25_tertile, dat2$death_in_hosp)

table(dat2$surgery_category, dat2$pm25_tertile)

dat2 %>%
  group_by(pm25_tertile) %>%
  summarise(mean_age = mean(age_years, na.rm=TRUE))

prop.trend.test(
  x = tapply(dat2$death_in_hosp, dat2$pm25_tertile, sum),
  n = tapply(dat2$death_in_hosp, dat2$pm25_tertile, length)
)

risk_tab <- dat2 %>%
  group_by(pm25_tertile) %>%
  summarise(
    n = n(),
    deaths = sum(death_in_hosp),
    risk = deaths / n
  )

risk_tab

risk_high <- risk_tab$risk[3]
risk_low  <- risk_tab$risk[1]
risk_high - risk_low

dat3 <- dat2 %>%
  left_join(imv_hours_total, by="hospitalization_id") %>%
  left_join(imv_hours_72, by="hospitalization_id")

kruskal.test(imv_hours_total ~ factor(pm25_tertile), data=dat3)
kruskal.test(imv_hours_total ~ factor(no2_tertile), data=dat3)

ggplot(dat3, aes(factor(pm25_tertile), imv_hours_total)) +
  geom_boxplot() +
  geom_jitter(width=0.1)

ggplot(dat3, aes(factor(no2_tertile), imv_hours_total)) +
  geom_boxplot() +
  geom_jitter(width=0.1)






