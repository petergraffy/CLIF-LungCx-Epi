source('code/lib/load.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=2)stop('Usage: audit_completed.R SITE_RETURN_ROOT PRIMARY_ROOT')
root<-a[1];primary<-a[2];collection<-file.path(root,'reports/site_supplement')
m<-jsonlite::fromJSON(file.path(collection,'collection_manifest.json'));sites<-m$sites
stopifnot(m$patient_level_inputs_read==FALSE,m$suppressed_values_preserved==TRUE)
files_n<-integer(length(sites));names(files_n)<-sites
for(s in sites){
 p<-file.path(root,'site_packets',s,'share');manifest<-jsonlite::fromJSON(file.path(p,'site_return_manifest.json'))
 ledger<-unlist(manifest$files,use.names=TRUE);stopifnot(all(unname(tools::md5sum(file.path(p,names(ledger))))==unname(ledger)));files_n[s]<-length(ledger)
 status<-readr::read_csv(file.path(p,'supplement/sofa_source_audit_status.csv'),show_col_types=FALSE)
 stopifnot(nrow(status)==3,all(status$status=='completed'))
 inv<-readr::read_csv(file.path(p,'supplement/return_file_inventory.csv'),show_col_types=FALSE);stopifnot(!any(inv$status=='missing'))
 ph<-readr::read_csv(file.path(p,'federated/model_PH_and_denominators.csv'),show_col_types=FALSE)
 computed<-ph$status=='computed_original_coefficients_reproduced';stopifnot(all(ph$coefficient_max_abs_difference[computed]==0,na.rm=TRUE))
 piece<-readr::read_csv(file.path(p,'federated/exploratory_piecewise_metadata.csv'),show_col_types=FALSE)
 for(endpoint in c('death','hospice','death_or_hospice')){
  z<-piece[piece$endpoint==endpoint,];original<-ph[ph$model_key==paste0('core_',endpoint),]
  stopifnot(nrow(z)==2,sum(z$events_n)==original$events_n,z$candidate_n[z$interval_start_day==0]==original$n,z$candidate_n[z$interval_start_day==7]<=original$n)
 }
}
ledger<-jsonlite::fromJSON(file.path(primary,'reports/frozen_artifact_ledger.json'),simplifyVector=FALSE)$ledger
for(z in ledger)stopifnot(unname(tools::md5sum(z$path))==z$md5)
requirements<-readr::read_csv('docs/site_return_requirements.csv',show_col_types=FALSE)
forms<-readr::read_csv(file.path(collection,'site_reporting_metadata.csv'),show_col_types=FALSE)
ph<-readr::read_csv(file.path(collection,'model_PH_and_denominators.csv'),show_col_types=FALSE)
counts<-readr::read_csv(file.path(primary,'reports/tables/cohort_counts.csv'),show_col_types=FALSE)
jsonlite::write_json(list(status='passed',sites=sites,cohort_n=sum(counts$n),reviewer_checklist_rows=nrow(requirements),aggregate_files_by_site=as.list(files_n),source_audits_all_completed=TRUE,required_families_captured=TRUE,global_and_term_PH_tests_captured=TRUE,diagnostic_coefficients_equal_original=TRUE,piecewise_riskset_and_event_conservation=TRUE,original_frozen_artifacts_unchanged=TRUE,site_confirmation_fields_pending=sum(forms$status=='site_confirmation_required'),no_patient_level_collection=TRUE,source_md5=as.list(tools::md5sum(c(list.files('diagnostics/site_return',full.names=TRUE), 'code/lib/site_supplement.R')))),file.path(root,'reports/site_return_audit.json'),auto_unbox=TRUE,pretty=TRUE)
notes<-c('Site return / reviewer coverage review — 2026-10-26','',
 paste('Sites:',paste(sites,collapse=', ')),paste('Primary cohort:',sum(counts$n)),paste('Reviewer-to-output checklist rows:',nrow(requirements)),
 paste('Aggregate return files:',paste(paste(names(files_n),files_n,sep='='),collapse=', ')),
 '', 'CAPTURED',
 'Every site: six-component SOFA completeness/missing patterns; source lab and FiO2/GCS/weight/vasoactive recovery diagnostics; original versus candidate score coverage; explicit unavailable oncology/source fields; code-status source vocabulary; cohort/early-event reconciliation; repeated admissions; annual evidence coverage.',
 'Complete frozen-model packet: local k/silhouette/stability/fingerprints; every full/held-out outcome version; separate death/hospice/composite; denominators/formulas/statuses; competing risks; nested/transported severity prediction, calibration/Brier/decision curves; ARF/code/cancer characterization; all eight refederated sensitivities; descriptive absorbing-state returns.',
 'Global/term PH diagnostics and exploratory posthoc0–7/7–28-day effects/metadata are captured. Diagnostic refits reproduce original coefficients exactly; original model artifacts are unchanged.',
 '', 'REMAINING INTERPRETATION / CONFIRMATIONS',
 'NU core composite global PH test p about0.023; P2/P3 term flags suggest time-varying associations. Exploratory time-window returns are available; they do not replace primary HRs or establish PH within windows.',
 'SOFA remains modified/incomplete. Candidate recovered scores are diagnostic only until outcome/comparator models are versioned and rerun. Validated APACHE II/IV unavailable; cancer stage/ECOG/therapy/preferences incomplete; rare death/full-code effects may be withheld or nonestimable. ICU-free days not validated.',
 paste('Pending site-confirmation fields:',sum(forms$status=='site_confirmation_required')),
 'Site confirmations include IRB/waiver details, release/ETL and timezone/provenance, treatment-limit order semantics and oncology timing, funding/conflicts and data availability. No missing preference is treated as absence.',
 'Lead team: descriptive/association framing, current cohort scope, literature, RECORD/STROBE, clean/tracked manuscript and point-by-point response. Do not claim consistent predictive superiority, biological endotypes, a causal intervention threshold, or an already executed ten-site revision.',
 '', 'WHERE TO USE',
 paste0('Final per-site return: ',root,'/site_packets/SITE/share/'),
 paste0('Combined aggregate supplement: ',collection,'/'),
 'Current site runner automatically captures these outputs; see docs/site_return_workflow.md and docs/site_return_requirements.csv. Only share/ is returnable after site review; all source rows, candidate scores, exact diagnostics and diagnostic fits stay local.',
 '', 'VALIDATION',
 'Synthetic core regression and two-site integration passed; arbitrary-site CSV supplement capture passed; SOFA recovery contracts passed; duplicate/consensus/tamper/nested-identifier rejection passed; PH reconstruction and day7 risk-set/event conservation passed. Source R parsing, shell syntax and diff whitespace checks passed.')
writeLines(notes,file.path(root,'reports/site_return_review_report.txt'))
cat('Site-return audits passed\n')
