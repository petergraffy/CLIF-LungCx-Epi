source('diagnostics/revision20261022/score_lib.R');source('code/lib/central.R')
root<-rev_root();out<-file.path(root,'reports','tables');dir.create(out,recursive=TRUE,showWarnings=FALSE)
expected<-c('phenotype_characteristics.csv','code_transitions.csv','outcomes_by_code.csv','arf_type_profiles.csv','clinical_table.csv','sensitivity_associations.csv','sensitivity_metadata.csv','early_events.csv','cohort_flow.csv','severity_availability.csv','arf_associations.csv','arf_joint_tests.csv','review_manifest.json')
for(site in fed_sites()){
 from<-file.path(root,'sites',site,'export');to<-file.path(root,'exchange','review',site);dir.create(to,recursive=TRUE,showWarnings=FALSE)
 for(nm in expected)stopifnot(file.copy(file.path(from,nm),file.path(to,nm),overwrite=identical(Sys.getenv('CLIF_REVISION_REFRESH_REVIEW'),'1')))
 for(nm in c('transport_metrics.csv','paired_increments.csv','decision_curves.csv','transport_manifest.json'))stopifnot(file.copy(file.path(from,'scores',nm),file.path(root,'exchange','scores',site,nm),overwrite=identical(Sys.getenv('CLIF_REVISION_REFRESH_REVIEW'),'1')))
}
for(nm in expected[grepl('csv$',expected)]){
 z<-dplyr::bind_rows(lapply(fed_sites(),function(s)readr::read_csv(file.path(root,'exchange','review',s,nm),show_col_types=FALSE)))
 fed_write(z,file.path(out,nm))
 if(nm=='sensitivity_associations.csv')fed_write(z |> dplyr::group_by(subset,model,endpoint,phenotype,reference) |> dplyr::group_modify(~pool_effect(.x)) |> dplyr::ungroup(),file.path(out,'pooled_adjusted_sensitivities.csv'))
 if(nm=='arf_associations.csv')fed_write(z |> dplyr::group_by(endpoint,term,effect) |> dplyr::group_modify(~pool_effect(.x)) |> dplyr::ungroup(),file.path(out,'pooled_arf_associations.csv'))
}
for(nm in c('transport_metrics.csv','paired_increments.csv','decision_curves.csv','model_metadata.csv')){
 z<-dplyr::bind_rows(lapply(fed_sites(),function(s)readr::read_csv(file.path(root,'exchange','scores',s,nm),show_col_types=FALSE)))
 fed_write(z,file.path(out,paste0('score_',nm)))
}
# Explicit availability row prevents a missing score being mistaken for a low AUC.
fed_write(tidyr::expand_grid(site=fed_sites(),score=c('APACHE II','APACHE IV')) |> dplyr::mutate(status='not_computable_from_current_extract',auc=NA_real_,reason='No recorded score; required chronic-health/admission context unverified; additional physiology unavailable. No APACHE proxy labeled as a validated score.'),file.path(out,'apache_comparison_status.csv'))
