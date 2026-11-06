# Validate the complete separate fallback run and its I/O/collection contract.
source('python_fallback/environment.R');clif_fallback_check_environment()
source('scripts/deployment_lib.R');source('python_fallback/verify_return.R')
a<-commandArgs(trailingOnly=TRUE);stopifnot(length(a)==1L)
root<-normalizePath(a[1]);workspace<-file.path(root,'local/pipeline');ret<-file.path(root,'analysis/stage2_return')
m<-deploy_validate_inventory(ret);b<-deploy_validate_inventory('frozen/stage2_v1','bundle_manifest.json');bhash<-unname(tools::md5sum('frozen/stage2_v1/bundle_manifest.json'))
verify_stage2_return(m,b,bhash)
stopifnot(!m$clustering_performed,!m$stage1_outputs_required,m$input_mode=='rebuild_from_clif',m$primary_consensus_id==b$primary_consensus_id)
for(f in c('profiles/hourly_joint.csv','profiles/initial_support_to_phenotype.csv','profiles/phenotype_transport_agreement.csv','prediction/score_histograms.csv','prediction/fold_audit.csv','survival/exit_histograms.csv','descriptives/categorical.csv','four_endpoint_cox/primary/associations.csv','supplement_package/supplement/manifest.json'))stopifnot(file.exists(file.path(ret,f)))
stopifnot(length(list.files(file.path(ret,'supplement_package'),recursive=TRUE))>20)
session<-readLines(file.path(ret,'environment/R_sessionInfo.txt'))
stopifnot(!any(grepl('(^|[[:space:]])arrow_[0-9]',session)),!'arrow'%in%loadedNamespaces())
for(variant in unlist(b$variants)){
 original<-file.path('frozen/stage2_v1',variant,'coordinator/full/prototypes.csv')
 copied<-file.path(workspace,'frozen/stage2_python_v1',variant,'coordinator/full/prototypes.csv')
 stopifnot(unname(tools::md5sum(original))==unname(tools::md5sum(copied)))
}
# Validate a real native packet alongside the fallback contract without changing either.
native<-'output/run_NU_stage2_sourceonly_20261105_v1_stage2_stage2_v1_independent/stage2_return'
if(dir.exists(native))verify_stage2_return(deploy_validate_inventory(native),b,bhash)
bad<-m;bad$source_hash<-'not_an_issued_source';stopifnot(inherits(tryCatch(verify_stage2_return(bad,b,bhash),error=identity),'error'))
bad<-m;bad$python_fallback$parent_bundle_md5<-'different_frozen_mapping';stopifnot(inherits(tryCatch(verify_stage2_return(bad,b,bhash),error=identity),'error'))
finder<-clif_package_matches
clif_package_matches<-function(package,version,...){if(package=='arrow')stop('R Arrow must not be checked');finder(package,version,...)}
clif_fallback_check_environment();clif_package_matches<-finder
# Exercise the shipped bridge, not the experimental prototype, without R Arrow.
python<-jsonlite::fromJSON(file.path(workspace,'config/config.json'))$python
Sys.setenv(CLIF_FALLBACK_PYTHON=python)
source('python_fallback/parquet.R')
old<-getwd();setwd(workspace)
fixture<-file.path(root,'edge_cases.parquet')
x<-data.frame(hospitalization_id=c('00001','9007199254740993','third','fourth'),text=c('',NA,'NA','comma, quote" newline\n雪'),value=c(pi,NA,NaN,-1e-100),flag=c(TRUE,FALSE,NA,TRUE),big=bit64::as.integer64(c('9007199254740993',NA,'-9007199254740993','42')),time=as.POSIXct('2020-01-01',tz='UTC')+c(0,.125,NA,3600),date=as.Date('2020-01-01')+c(0,1,NA,2))
clif_write_parquet(x,fixture);y<-clif_read_parquet(fixture)
stopifnot(isTRUE(all.equal(x,y,check.attributes=FALSE,tolerance=0)))
z<-clif_read_parquet(fixture,columns=c('hospitalization_id','value'),id_column='hospitalization_id',ids='9007199254740993')
stopifnot(nrow(z)==1L,z$hospitalization_id=='9007199254740993',is.na(z$value))
z<-clif_read_parquet(fixture,columns='hospitalization_id',id_column='hospitalization_id',ids=character());stopifnot(nrow(z)==0L)
distinct<-file.path(root,'distinct.parquet');clif_write_parquet(data.frame(category=c('full','full',NA,'dni')),distinct)
stopifnot(nrow(clif_read_parquet(distinct,columns='category',distinct=TRUE))==3L)
setwd(old);stopifnot(!'arrow'%in%loadedNamespaces())
cat('PASS: complete independent Stage2 return; all export families, templates, native/fallback verification, tamper rejection, ignored R Arrow version, and exact shipped Parquet transport without R Arrow.\n')
