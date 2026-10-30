source('diagnostics/revision20261022/lib.R')
root<-'output/run_three_site_refresh_20261026';dir.create(file.path(root,'reports/tables'),recursive=TRUE,showWarnings=FALSE)
for(site in c('UCMC','NU','RUSH')){
 status<-system2(file.path(R.home('bin'),'Rscript'),c('diagnostics/refresh20261026/site_compare.R',site),stdout=file.path(root,'logs',paste0('compare_',site,'.log')),stderr=TRUE);stopifnot(status==0)
}
for(nm in c('source_comparison.csv','assignment_transitions.csv','partition_comparison.csv')){
 rows<-list()
 for(site in c('UCMC','NU','RUSH')){
  from<-file.path(root,'sites',site,'export',nm);to<-file.path(root,'exchange',site,nm);dir.create(dirname(to),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(from,to,overwrite=FALSE));rows[[site]]<-readr::read_csv(to,show_col_types=FALSE)
 }
 fed_write(dplyr::bind_rows(rows),file.path(root,'reports/tables',nm))
}
# Comparing summaries reads central aggregate reports only.
for(kind in c('previous','refreshed')){
 p<-if(kind=='previous')'output/run_federated_rehearsal_20261022' else old_root()
 for(nm in c('pooled_associations.csv','transport_metrics.csv','transport_coverage.csv','cohort_counts.csv')){
  z<-readr::read_csv(file.path(p,'reports/tables',nm),show_col_types=FALSE);fed_write(z,file.path(root,'reports/tables',paste0(kind,'_',nm)))
 }
}
