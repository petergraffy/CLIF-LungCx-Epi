# End-to-end one-pass smoke from a previously generated synthetic Stage1 fixture.
source('scripts/deployment_lib.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=2L)stop('Usage: tests/stage2_site.R SYNTHETIC_STAGE1_CONFIG FROZEN_BUNDLE')
config<-normalizePath(a[1]);cfg<-read_run_config(config);original<-deploy_root(cfg)
stopifnot(cfg$site_name=='SYNTHETIC_SITE',grepl('run_deployment_site_',cfg$analysis$run_id))
base<-tempfile('run_stage2_site_',tmpdir='output');dir.create(base);base<-normalizePath(base)
bundle<-file.path(base,paste0('synthetic_',basename(base)));log<-file.path(base,'admit.log')
s<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','scripts/admit_stage2_site.R',shQuote(a[2]),shQuote(file.path(original,'stage1_return')),shQuote(bundle),shQuote('Synthetic fixture admission for a one-pass integration test.')),stdout=log,stderr=log)
if(s!=0){cat(readLines(log),sep='\n');stop('admission failed')}
# Unknown sites must still fail before any copying or modeling.
bm<-deploy_validate_inventory(bundle,'bundle_manifest.json');bad<-bm;bad$sites<-list('OTHER');stopifnot(inherits(tryCatch(deploy_stage2_context(cfg,bundle,bad,FALSE),error=identity),'error'))
source_before<-unname(tools::md5sum(file.path(original,'stage1_return/return_manifest.json')))
for(command in c('stage2-preflight','stage2')){
 log<-file.path(base,paste0(command,'.log'));s<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','scripts/site.R',command,shQuote(config),shQuote(bundle)),stdout=log,stderr=log)
 if(s!=0){cat(tail(readLines(log),25),sep='\n');stop('One-pass smoke failed; inspect ',base)}
}
context<-deploy_stage2_context(cfg,bundle,bm,FALSE);ret<-file.path(context$root,'stage2_return');m<-deploy_validate_inventory(ret)
stopifnot(m$stage=='stage2',m$status=='complete',!m$external_site_outcome_validation,m$validation=='patient_grouped_local_cross_validation',source_before==unname(tools::md5sum(file.path(original,'stage1_return/return_manifest.json'))))
for(f in c('prediction/score_histograms.csv','prediction/fold_audit.csv','survival/exit_histograms.csv','profiles/hourly_joint.csv','profiles/initial_support_to_phenotype.csv','profiles/phenotype_transport_agreement.csv','descriptives/categorical.csv','four_endpoint_cox/primary/associations.csv'))stopifnot(file.exists(file.path(ret,f)))
f<-readr::read_csv(file.path(ret,'prediction/fold_audit.csv'),show_col_types=FALSE);stopifnot(all(f$patient_overlap==0))
# Exercise exact central survival reconstruction without another site pass.
log<-file.path(base,'central.log');s<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','diagnostics/stage2/central_survival.R',shQuote(file.path(base,'central_survival')),shQuote(file.path(ret,'survival'))),stdout=log,stderr=log)
if(s!=0){cat(tail(readLines(log),20),sep='\n');stop('central histogram reconstruction failed')}
s<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','diagnostics/pooled_survival/central_pool.R',shQuote(file.path(base,'pooled_survival')),shQuote(file.path(base,'central_survival/grid')),shQuote(file.path(base,'central_survival/sites/SYNTHETIC_SITE'))),stdout=log,stderr=log)
if(s!=0){cat(tail(readLines(log),20),sep='\n');stop('central AJ/KM pooling failed')}
s<-system2(cfg$python,c(shQuote('diagnostics/stage2/central_prediction.py'),shQuote(file.path(base,'pooled_prediction')),shQuote(file.path(ret,'prediction'))),stdout=log,stderr=log)
if(s!=0){cat(tail(readLines(log),20),sep='\n');stop('central ROC/tdAUC pooling failed')}
cat('PASS: late-site admission preserves full prototypes; one site command exports all models/sensitivities/descriptors/curves/prediction; central AJ/KM/ROC/tdAUC from aggregates only.\nReturn: ',ret,'\nTest logs: ',base,'\n',sep='')
