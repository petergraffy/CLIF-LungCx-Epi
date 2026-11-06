# End-to-end Stage2 smoke from source CLIF fixture tables only; no Stage1 outputs.
source('scripts/deployment_lib.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=2L)stop('Usage: tests/stage2_site.R SYNTHETIC_STAGE1_CONFIG FROZEN_BUNDLE')
fixture_config<-normalizePath(a[1]);cfg<-read_run_config(fixture_config)
stopifnot(cfg$site_name=='SYNTHETIC_SITE',grepl('run_deployment_site_',cfg$analysis$run_id))
base<-tempfile('run_stage2_site_',tmpdir='output');dir.create(base);base<-normalizePath(base)
cfg$analysis$run_id<-basename(base);cfg$stage1_root<-file.path(base,'DELETED_STAGE1');config<-file.path(base,'site.json');deploy_json(cfg,config)
stopifnot(!dir.exists(deploy_root(cfg)),!dir.exists(init_path<-file.path(cfg$repo,'output',paste0('run_',cfg$site_name,'_',cfg$analysis$run_id))))
bundle<-file.path(base,'independent_bundle');deploy_tree(a[2],bundle,c('csv','json','md'));bm<-deploy_read(file.path(bundle,'bundle_manifest.json'))
for(v in unlist(bm$variants)){
 dest<-file.path(bundle,v,'coordinator/loo_SYNTHETIC_SITE');deploy_tree(file.path(bundle,v,'coordinator/full'),dest,c('csv','json'))
 fm<-deploy_read(file.path(dest,'frozen_manifest.json'));fm$version<-'loo_SYNTHETIC_SITE';fm$held_out_site<-'SYNTHETIC_SITE';deploy_json(fm,file.path(dest,'frozen_manifest.json'))
}
bm$sites<-c(bm$sites,list('SYNTHETIC_SITE'));bm$workspace_id<-basename(base);bm$files<-deploy_inventory(bundle);deploy_json(bm,file.path(bundle,'bundle_manifest.json'))
protos<-list.files(bundle,pattern='prototypes.csv',recursive=TRUE,full.names=TRUE);before<-tools::md5sum(protos)
bad<-bm;bad$sites<-list('OTHER');stopifnot(inherits(tryCatch(deploy_stage2_context(cfg,bundle,bad,FALSE),error=identity),'error'))
for(command in c('stage2-preflight','stage2')){
 log<-file.path(base,paste0(command,'.log'));s<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','scripts/site.R',command,shQuote(config),shQuote(bundle)),stdout=log,stderr=log)
 if(s!=0){cat(tail(readLines(log),25),sep='\n');stop('One-pass smoke failed; inspect ',base)}
}
context<-deploy_stage2_context(cfg,bundle,bm,FALSE);ret<-file.path(context$root,'stage2_return');m<-deploy_validate_inventory(ret)
stopifnot(m$stage=='stage2',m$status=='complete',!m$external_site_outcome_validation,m$validation=='patient_grouped_local_cross_validation',identical(before,tools::md5sum(protos)),!dir.exists(deploy_root(cfg)),!dir.exists(init_path),!dir.exists(cfg$stage1_root))
for(f in c('input_audit/cohort_rebuild_comparison.csv','prediction/score_histograms.csv','prediction/fold_audit.csv','survival/exit_histograms.csv','profiles/hourly_joint.csv','profiles/initial_support_to_phenotype.csv','profiles/phenotype_transport_agreement.csv','descriptives/categorical.csv','four_endpoint_cox/primary/associations.csv'))stopifnot(file.exists(file.path(ret,f)))
f<-readr::read_csv(file.path(ret,'prediction/fold_audit.csv'),show_col_types=FALSE);stopifnot(all(f$patient_overlap==0))
# Exercise exact central survival reconstruction without another site pass.
log<-file.path(base,'central.log');s<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','diagnostics/stage2/central_survival.R',shQuote(file.path(base,'central_survival')),shQuote(file.path(ret,'survival'))),stdout=log,stderr=log)
if(s!=0){cat(tail(readLines(log),20),sep='\n');stop('central histogram reconstruction failed')}
s<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','diagnostics/pooled_survival/central_pool.R',shQuote(file.path(base,'pooled_survival')),shQuote(file.path(base,'central_survival/grid')),shQuote(file.path(base,'central_survival/sites/SYNTHETIC_SITE'))),stdout=log,stderr=log)
if(s!=0){cat(tail(readLines(log),20),sep='\n');stop('central AJ/KM pooling failed')}
s<-system2(cfg$python,c(shQuote('diagnostics/stage2/central_prediction.py'),shQuote(file.path(base,'pooled_prediction')),shQuote(file.path(ret,'prediction'))),stdout=log,stderr=log)
if(s!=0){cat(tail(readLines(log),20),sep='\n');stop('central ROC/tdAUC pooling failed')}
cat('PASS: no Stage1 folders exist; all inputs rebuilt with clustering disabled and frozen templates unchanged; one site command exports all models/sensitivities/descriptors/curves/prediction; central AJ/KM/ROC/tdAUC from aggregates only.\nReturn: ',ret,'\nTest logs: ',base,'\n',sep='')
