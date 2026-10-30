source('code/lib/load.R')
args<-commandArgs(trailingOnly=TRUE);if(length(args)!=3)stop('Usage: model_diagnostics.R SITE PRIMARY_ROOT NEW_OUTPUT')
site<-args[1];home<-file.path(args[2],'sites',site);out<-args[3]
if(dir.exists(out))stop('Choose fresh diagnostic output')
dir.create(file.path(out,'local'),recursive=TRUE);dir.create(file.path(out,'share'),recursive=TRUE)
saved<-readRDS(file.path(home,'local/full/outcomes.rds'));metadata<-readr::read_csv(file.path(home,'export/outcomes/full/model_metadata.csv'),show_col_types=FALSE)
z<-supplement_model_diagnostics(saved,metadata);saveRDS(z,file.path(out,'local/diagnostic_fits.rds'));readr::write_csv(z$table,file.path(out,'local/exact.csv'))
x<-protect_aggregate(z$table,c('n','events_n','non_events_n'),character(),10)
x$site<-site;readr::write_csv(x,file.path(out,'share/model_PH_and_denominators.csv'),na='')
terms<-protect_aggregate(z$term_tests,c('n','events_n','non_events_n'),'model_key',10);terms$site<-site;readr::write_csv(terms,file.path(out,'share/model_PH_terms.csv'),na='')
stopifnot(!any(z$table$status%in%c('diagnostic_reconstruction_mismatch','diagnostic_refit_failed','PH_test_unavailable')))
jsonlite::write_json(list(site=site,status='passed',original_coefficients_reproduced=TRUE,original_model_artifacts_unchanged=TRUE,diagnostic_refit_only=TRUE,nonestimable_original_models_not_tested=TRUE,source_md5=as.list(tools::md5sum(c('diagnostics/site_return/model_diagnostics.R','code/lib/site_supplement.R')))),file.path(out,'share/manifest.json'),auto_unbox=TRUE,pretty=TRUE)
piece<-supplement_piecewise_associations(saved$data,10)
saveRDS(piece,file.path(out,'local/piecewise_diagnostic_fits.rds'))
a<-piece$associations;a$site<-site
readr::write_csv(a,file.path(out,'share/exploratory_piecewise_associations.csv'),na='')
b<-protect_aggregate(piece$metadata,c('candidate_n','fitted_n','events_n','non_events_n'),character(),10);b$site<-site
readr::write_csv(b,file.path(out,'share/exploratory_piecewise_metadata.csv'),na='')
print(z$table)
