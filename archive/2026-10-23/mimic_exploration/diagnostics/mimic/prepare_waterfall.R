source('diagnostics/mimic/portable_lib.R')
root<-mimic_root();a<-readRDS(file.path(root,'local/prepared.rds'));d<-a$landmark[a$landmark$landmark_reason=='eligible',]
stopifnot(all(d$time_from_landmark>0))
b<-d |> dplyr::transmute(job='h72_primary',hospitalization_id,trajectory_origin=t0,start=pmax(admission_dttm,t0-6*3600),end=t0+72*3600)
rs<-a$tables$respiratory_support;rs$fio2_set<-standard_fio2(rs$fio2_set)
arrow::write_parquet(rs,file.path(root,'local/raw_respiratory.parquet'));arrow::write_parquet(b,file.path(root,'local/window_bounds.parquet'))
fed_json(list(site='MIMIC',boundary_rule='six-hour same-admission lookback; input/output strictly before ICU+72h; alive/in-hospital landmark cohort',numeric_bfill=FALSE,original_physiology_preserved=TRUE),file.path(root,'local/input_manifest.json'))
