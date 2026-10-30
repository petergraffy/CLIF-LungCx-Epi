source('code/lib/load.R');source('diagnostics/manuscript_figures/survival_tdauc_lib.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=3)stop('Usage: site_survival_tdauc.R SITE PRIMARY_ROOT NEW_OUTPUT')
site<-a[1];primary<-a[2];out<-a[3];if(dir.exists(out))stop('Choose a fresh site output')
home<-file.path(primary,'sites',site);disc<-readRDS(file.path(home,'local/discovery.rds'))
input<-jsonlite::fromJSON(file.path(disc$input_root,'local/input_manifest.json'));prepared<-readRDS(file.path(input$baseline_root,'local/prepared.rds'))
all<-prepared$cohort
# All eligible admissions have observable discharge timestamps in these extracts.
time<-hours_between(all$event_time,all$t0)/24
stopifnot(all(is.finite(time)&time>0))
km<-km_release(all,time) |> dplyr::mutate(clock='icu_admission',site=site)
full<-readRDS(file.path(home,'local/full/outcomes.rds'))$data
for(ph in sort(unique(full$phenotype))){z<-full[full$phenotype==ph,];km<-dplyr::bind_rows(km,km_release(z,z$time_from_landmark,ph) |> dplyr::mutate(clock='landmark72h',site=site))}
dir.create(file.path(out,'local'),recursive=TRUE);share<-file.path(out,'share');dir.create(share)
readr::write_csv(km,file.path(share,'survival_curves.csv'),na='')
audit<-tibble::tibble(site=site,cohort=c('eligible','landmark'),n=c(nrow(all),nrow(full)),unknown_outcome_n=c(sum(!all$valid_outcome),sum(!full$valid_outcome)))
audit<-protect_aggregate(audit,c('n','unknown_outcome_n'),'cohort',10L)
readr::write_csv(audit,file.path(share,'survival_denominators.csv'),na='')
loo<-readRDS(file.path(home,'local',paste0('loo_',site),'outcomes.rds'))$data
pr<-readRDS(file.path(home,'local',paste0('loo_',site),'transport_extended_private.rds'))
stopifnot(identical(pr$patient_id,loo$patient_id))
metrics<-list();roc<-list()
capture<-function(d,p,cohort,endpoint){
 x<-td_capture(d,p,cohort,endpoint,site)
 metrics[[paste(cohort,endpoint)]]<<-x$metrics;roc[[paste(cohort,endpoint)]]<<-x$roc
 cat(site,cohort,endpoint,'tdAUC complete\n');flush.console()
}
capture(loo,pr$predictions,'full_landmark','death_or_hospice')
aproot<-Sys.getenv('CLIF_APACHE_OUTPUT',unset='output/run_APACHE_adapted_20261030_v2')
ap<-readRDS(file.path(aproot,'sites',site,'local/evaluation.rds'))
revroot<-Sys.getenv('CLIF_REVISION_OUTPUT',unset='output/run_revision_gaps_20261026')
so<-readRDS(file.path(revroot,'sites',site,'local/score_transport.rds'))
for(endpoint in c('death','death_or_hospice')){
 # Full-cohort death-specific frozen predictors already exist in the AP extension.
 if(endpoint=='death'){
  x<-ap[['all_case_missing_indicators death']];d<-loo[match(x$ids,loo$hospitalization_id),];stopifnot(identical(d$hospitalization_id,x$ids))
  capture(d,x$pred[c('baseline','phenotype','summaries','summaries_phenotype')],'full_landmark',endpoint)
 }
 x<-so[[paste('measured_sofa',endpoint)]];d<-loo[match(x$ids,loo$hospitalization_id),];stopifnot(identical(d$hospitalization_id,x$ids))
 capture(d,x$predictions[c('baseline','phenotype','sofa','sofa_summaries','sofa_summaries_phenotype')],'measured_sofa',endpoint)
 for(cohort in c('complete_hb_proxy','matched_complete_SOFA')){
  x<-ap[[paste(cohort,endpoint)]];d<-loo[match(x$ids,loo$hospitalization_id),];stopifnot(identical(d$hospitalization_id,x$ids))
  kinds<-if(cohort=='complete_hb_proxy')c('baseline','phenotype','baseline_adapted','adapted_summaries','adapted_summaries_phenotype') else c('phenotype','baseline_sofa','baseline_adapted','sofa_summaries_phenotype','adapted_summaries_phenotype')
  capture(d,x$pred[kinds],cohort,endpoint)
 }
}
metric<-dplyr::bind_rows(metrics);rocs<-dplyr::bind_rows(roc)
# Primary AUC_2 at28d must reproduce each frozen binary-horizon AUROC.
old<-readr::read_csv(file.path(primary,'reports/tables/transport_metrics.csv'),show_col_types=FALSE) |> dplyr::filter(site==.env$site,metric=='auc',model%in%names(pr$predictions))
new<-metric |> dplyr::filter(cohort=='full_landmark',endpoint=='death_or_hospice',day==28,control_definition==2)
check<-dplyr::inner_join(old,new,by='model',suffix=c('_old','_new'));stopifnot(nrow(check)==length(pr$predictions),all(abs(check$estimate_old-check$estimate_new)<1e-10))
for(pair in list(list(x=metric,name='tdauc.csv'),list(x=rocs,name='tdroc.csv'))){
 stopifnot(!any(names(pair$x)%in%c('patient_id','hospitalization_id','t0','event_time','marker','threshold')))
 readr::write_csv(pair$x,file.path(share,pair$name),na='')
}
files<-list.files(share,full.names=TRUE)
jsonlite::write_json(list(status='passed',site=site,primary_consensus_id=jsonlite::fromJSON(file.path(primary,'coordinator/full/frozen_manifest.json'))$consensus_id,
 aggregate_policy='exact counts; no disclosure suppression',whole_cohort_clock='ICU admission; all eligible stays; unknown disposition censored at observed discharge with explicit missing-outcome assumption',
 phenotype_clock='72h landmark; no retrospective baseline phenotype assignment',
 survival='Cause-specific KM; competing hospital exits censored; not all-cause survival after discharge; CIF retained separately',
 timeROC_version=as.character(utils::packageVersion('timeROC')),tdauc='Cumulative/dynamic timeROC IPCW; control2 includes competing events; control1 only event-free controls',
 no_ordinary_censoring_before28horizon=TRUE,IPCW_unit_weights_verified=TRUE,bootstrap='200 patient-cluster resamples; pointwise intervals conditional on frozen predictors',
 prediction='Fixed frozen28day risk markers evaluated over time; no horizon-specific refits',day28_primary_auc_reproduced=TRUE,
 tdroc='Fixed FPR grid at7/14/28d; no score cutoffs exported',clinical_model_refits=FALSE,
 files=as.list(setNames(unname(tools::md5sum(files)),basename(files)))),file.path(share,'manifest.json'),auto_unbox=TRUE,pretty=TRUE)
cat(site,'survival and tdAUC audit passed\n')
