source('code/lib/load.R')
out<-Sys.getenv('CLIF_UCMC_REFRESH_OUTPUT',unset='output/run_UCMC_etl_refresh_audit_20261026');a<-readRDS(file.path(out,'local/paired_audit.rds'));path<-a$cfg$tables_path
columns<-list(patient=c('patient_id','sex_category','race_category','ethnicity_category','death_dttm'),hospitalization=c('patient_id','hospitalization_id','admission_dttm','discharge_dttm','age_at_admission','discharge_category'),adt=c('hospitalization_id','in_dttm','out_dttm','location_category'),hospital_diagnosis=c('hospitalization_id','diagnosis_code','diagnosis_code_format','poa_present'))
t<-lapply(names(columns),function(n){x<-arrow::open_dataset(file.path(path,paste0('clif_',n,'.parquet'))) |> dplyr::select(dplyr::all_of(columns[[n]])) |> dplyr::collect();for(ts in grep('_dttm$',names(x),value=TRUE))x[[ts]]<-utc(x[[ts]]);x});names(t)<-names(columns)
new<-build_cohort(t,a$cfg)$cohort;old<-a$data;new<-new[match(old$hospitalization_id,new$hospitalization_id),]
stopifnot(!anyNA(new$hospitalization_id));vars<-c('t0','admission_dttm','discharge_dttm','age_years','event_type','event_time','valid_outcome','sex_category','race_category','ethnicity_category')
rows<-lapply(vars,function(v){x<-old[[v]];y<-new[[v]];different<-is.na(x)!=is.na(y) | (!is.na(x)&!is.na(y)&x!=y);tibble::tibble(variable=v,changed_hospitalizations=sum(different,na.rm=TRUE))})
readr::write_csv(dplyr::bind_rows(rows),file.path(out,'aggregate/anchor_field_changes.csv'));saveRDS(new,file.path(out,'local/current_matched_cohort.rds'))
