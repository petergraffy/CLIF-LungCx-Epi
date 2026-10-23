# Regression checks for date-only boundaries and optional code-status IDs.
source("code/lib/load.R")
check <- function(ok, label) {stopifnot(isTRUE(ok));cat("PASS:",label,"\n")}
run_input_checks <- function() {
 check(identical(utc("2018-01-01"),utc("2018-01-01 00:00:00")),"date-only study boundary parses at UTC midnight")
 mixed<-utc(c("2018-01-01",NA,"2020-06-02 12:34:56"))
 check(!is.na(mixed[1]) && is.na(mixed[2]) && format(mixed[3],tz="UTC")=="2020-06-02 12:34:56","mixed timestamps retain missingness and full times")
 t<-as.POSIXct(c("2017-12-31 12:00:00","2020-01-02 00:00:00","2025-01-01 00:00:00"),tz="UTC")
 ids<-c("before","within","after")
 tables<-list(patient=tibble::tibble(patient_id=paste0("p",1:3),sex_category="Female",race_category="White"),
 hospitalization=tibble::tibble(patient_id=paste0("p",1:3),hospitalization_id=ids,admission_dttm=t,discharge_dttm=t+10*86400,discharge_category="Home",age_at_admission=65),
 adt=tibble::tibble(hospitalization_id=ids,in_dttm=t+3600,out_dttm=t+4*86400,location_category="icu"),
 hospital_diagnosis=tibble::tibble(hospitalization_id=ids,diagnosis_code="C34.90",diagnosis_code_format="ICD10CM",poa_present=1))
 cfg<-list(site_name="SyntheticInput",repo=getwd(),analysis=run_defaults(),file_type="parquet")
 check(identical(build_cohort(tables,cfg)$cohort$hospitalization_id,"within"),"date-only study window excludes pre-window and end-boundary admissions")
 tables$respiratory_support<-tibble::tibble(hospitalization_id="within",recorded_dttm=t[2]+3600,device_category="Room Air")
 tables$vitals<-tibble::tibble(hospitalization_id="within",recorded_dttm=t[2]+3600,vital_category="spo2",vital_value=98)
 tables$labs<-tibble::tibble(hospitalization_id="within",lab_result_dttm=t[2]+3600,lab_category="creatinine",lab_value_numeric=.8)
 tables$medication_admin_continuous<-tibble::tibble(hospitalization_id="within",admin_dttm=t[2]+3600,med_category="norepinephrine",med_dose=0,med_dose_unit="mcg/kg/min")
 tables$patient_assessments<-tibble::tibble(hospitalization_id="within",recorded_dttm=t[2]+3600,assessment_category="gcs_total",numerical_value=15)
 tables$code_status<-tibble::tibble(hospitalization_id=c("within","foreign"),patient_id=c("p2","p999"),start_dttm=t[2]+c(0,1),code_status_category=c("Full","DNR"))
 cfg$tables_path<-tempfile("clif-input-check-");dir.create(cfg$tables_path)
 on.exit(unlink(cfg$tables_path,recursive=TRUE),add=TRUE)
 for(nm in names(tables))arrow::write_parquet(tables[[nm]],file.path(cfg$tables_path,paste0("clif_",nm,".parquet")))
 config_file<-file.path(cfg$tables_path,"test_config.json")
 jsonlite::write_json(cfg,config_file,auto_unbox=TRUE)
 check(!is.null(read_run_config(config_file)),"valid date-only run configuration is accepted")
 bad<-cfg;bad$analysis$start_date<-"invalid-date"
 jsonlite::write_json(bad,config_file,auto_unbox=TRUE)
 check(inherits(tryCatch(read_run_config(config_file),error=identity),"error"),"invalid study boundary is rejected rather than silently disabling filtering")
 bad<-cfg;bad$analysis$end_date<-"2017-01-01"
 jsonlite::write_json(bad,config_file,auto_unbox=TRUE)
 check(inherits(tryCatch(read_run_config(config_file),error=identity),"error"),"reversed study boundaries are rejected")
 loaded<-load_clif(cfg)
 check(nrow(loaded$code_status)==1 && loaded$code_status$patient_id=="p2","extra code-status hospitalization ID does not break projected patient filtering")
 check(identical(sort(names(loaded$code_status)),sort(c("patient_id","start_dttm","code_status_category"))),"optional identifier does not change clinical code-status fields")
}
run_input_checks()
