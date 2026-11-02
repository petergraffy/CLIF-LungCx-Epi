# Full portable Stage 1 smoke test on synthetic CLIF parquet data.
source('scripts/deployment_lib.R')
python<-Sys.getenv('CLIF_TEST_WATERFALL_PYTHON')
if(!nzchar(python))stop('Set CLIF_TEST_WATERFALL_PYTHON to the pinned Python 3.12 executable.')
# Reuse the clinical test fixture without executing its regression suite.
for(expr in parse('tests/run_tests.R'))if(is.call(expr)&&identical(expr[[1]],as.name('<-'))&&identical(expr[[2]],as.name('fixture')))eval(expr)
stopifnot(exists('fixture'))
base<-tempfile('run_deployment_site_',tmpdir='output');dir.create(base);base<-normalizePath(base)
tables<-file.path(base,'tables');dir.create(tables)
x<-fixture(80)
x$respiratory_support$device_category<-c('Room Air'='room air','IMV'='imv','High Flow NC'='high flow nc','Nasal Cannula'='nasal cannula')[x$respiratory_support$device_category]
x$respiratory_support$device_name<-x$respiratory_support$device_category
x$respiratory_support$lpm_set<-ifelse(x$respiratory_support$device_category=='room air',0,4)
for(n in names(x))arrow::write_parquet(x[[n]],file.path(tables,paste0('clif_',n,'.parquet')))
cfg<-jsonlite::fromJSON('config/config_template.json');cfg$site_name<-'SYNTHETIC_SITE';cfg$repo<-clif_root;cfg$tables_path<-tables;cfg$python<-if(grepl("^/|^[A-Za-z]:",python))python else file.path(clif_root,python);cfg$analysis$run_id<-basename(base)
config<-file.path(base,'site.json');deploy_json(cfg,config)
for(command in c('preflight','stage1')){
 log<-file.path(base,paste0(command,'.log'))
 status<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','scripts/site.R',command,shQuote(config)),stdout=log,stderr=log)
 if(status!=0){cat(tail(readLines(log),15),sep='\n');stop(command,' failed; diagnostics in ',base)}
}
root<-deploy_root(read_run_config(config));share<-file.path(root,'stage1_return');m<-deploy_validate_inventory(share)
stopifnot(m$status=='complete',m$stage=='stage1',!m$shared_phenotypes_assigned,!m$local_outcome_models_fit)
stopifnot(setequal(unlist(m$variants),DEPLOY_VARIANTS),file.exists(file.path(share,'descriptives/categorical.csv')))
stopifnot(!dir.exists(file.path(root,'primary/coordinator')),!dir.exists(file.path(root,'primary/sites/SYNTHETIC_SITE/export/outcomes')))
cat('PASS: full Stage 1 on synthetic source parquet; ten discoveries, descriptors, supplement, environment inventory and mandatory stop\n')
# An empty coordinator directory must not allow Stage 2.
bad<-file.path(base,'not_a_bundle');dir.create(bad)
status<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','scripts/site.R','stage2-fit',shQuote(config),shQuote(bad)),stdout=file.path(base,'bad_bundle.log'),stderr=file.path(base,'bad_bundle.log'))
stopifnot(status!=0)
cat('PASS: Stage 2 cannot proceed without a coordinator bundle\n')
cat('Synthetic Stage 1 results: ',share,'\n',sep='')
