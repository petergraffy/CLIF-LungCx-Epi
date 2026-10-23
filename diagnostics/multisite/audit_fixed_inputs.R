source("code/lib/load.R")
root<-"output/run_multisite_comparison_20261022"
manifest<-jsonlite::fromJSON(file.path(root,"local","comparison_manifest.json"),simplifyVector=FALSE)
audit<-list()
for(site in names(manifest$sites)) {
 baseline_root<-manifest$sites[[site]]$baseline_root
 current<-readRDS(file.path(baseline_root,"local","prepared.rds"))
 fitted<-readRDS(file.path(baseline_root,"local","derivation.rds"))$analyses$h72_primary
 # Rush's aborted preparation contains the extra 2017/2025 admissions and
 # allows verification that the fixed date-only parser reproduces v3c.
 source_prepared<-if(site=="RUSH")readRDS("output/run_RUSH_ccx_pilot_20261022_v3b_joint_primary/local/prepared.rds") else current
 cfg<-list(analysis=current$parameters);cfg$analysis$start_date<-"2018-01-01";cfg$analysis$end_date<-"2025-01-01"
 rebuilt<-build_cohort(source_prepared$tables,cfg)$cohort
 lm<-landmark_cohort(rebuilt,72)
 stopifnot(setequal(rebuilt$hospitalization_id,current$cohort$hospitalization_id),
 setequal(lm$hospitalization_id[lm$landmark_reason=="eligible"],fitted$data$hospitalization_id))
 audit[[site]]<-list(eligible_hospitalizations=nrow(rebuilt),primary_hospitalizations=nrow(fitted$data),
 fixed_date_only_parser_matches_completed_cohorts=TRUE,completed_analysis_source_id=manifest$sites[[site]]$baseline_source_id)
}
jsonlite::write_json(list(status="passed",sites=audit,fixed_core_source_id=code_fingerprint(),
 source_snapshot="analytical_source_before_input_fixes.tar.gz",input_regression_tests="input_contract_tests.log",synthetic_tests="synthetic_tests.log"),
 file.path(root,"local","input_fixes_audit.json"),auto_unbox=TRUE,pretty=TRUE)
message("Fixed date-only parsing reproduces all three completed cohorts and landmark sets")
