# Sourced only by scripts/site.R after environment, bundle and lineage validation.
if(!identical(bm$workflow,'one_pass_local_validation_v1'))stop('This bundle does not implement one-pass Stage2.')
for(v in DEPLOY_VARIANTS){
 dest<-file.path(deploy_analysis_root(root,v),'coordinator')
 deploy_tree(file.path(bundle,v,'coordinator'),dest,c('csv','json'))
 for(version in c('full',paste0('loo_',sites))){p<-deploy_load_prototypes(file.path(dest,version));saveRDS(p,file.path(dest,version,'prototypes.rds'))}
}
for(v in setdiff(DEPLOY_VARIANTS,'absorbing')){
 Sys.setenv(CLIF_FEDERATED_OUTPUT=deploy_analysis_root(root,v),CLIF_REVISION_VARIANT=v)
 script<-if(v=='primary')'diagnostics/federated/site_outcomes.R' else 'diagnostics/revision20261022/site_outcomes.R'
 # LOO evaluates phenotype transportability only; outcome validation is local CV.
 for(version in c('full',paste0('loo_',site)))step(paste('fit',v,version,sep='_'),script,c(site,version))
}
Sys.setenv(CLIF_FEDERATED_OUTPUT=file.path(root,'primary'))
for(script in c('site_review','site_descriptive_tables','site_arf_effects','site_severity_source_audit','site_code_source_audit'))step(script,paste0('diagnostics/revision20261022/',script,'.R'),site)
step('apache_prepare','diagnostics/apache_adapted/site_prepare.R',site)
step('local_validation','diagnostics/stage2/site_prediction.R',c(site,bundle,file.path(root,'prediction')))
step('survival','diagnostics/stage2/site_survival.R',c(site,bundle,file.path(root,'survival_share')))
Sys.setenv(CLIF_FEDERATED_OUTPUT=deploy_analysis_root(root,'absorbing'))
step('absorbing_summary','diagnostics/revision20261022/site_absorbing_summary.R',site)
step('absorbing_audit','diagnostics/revision20261022/site_absorbing_audit.R',site)
Sys.setenv(CLIF_FEDERATED_OUTPUT=file.path(root,'primary'))
step('four_cox','diagnostics/four_endpoint_cox/site_capture.R',c(site,file.path(root,'primary'),file.path(root,'revision'),file.path(root,'four_cox')))
step('complete_capture','diagnostics/stage2/site_capture.R',c(site,bundle,config,root))
return<-file.path(root,'stage2_return')
deploy_finish(return,'stage2',site,list(workflow=bm$workflow,primary_consensus_id=bm$primary_consensus_id,parameter_bundle_md5=unname(tools::md5sum(file.path(bundle,'bundle_manifest.json'))),validation='patient_grouped_local_cross_validation',external_site_outcome_validation=FALSE))
cat('STAGE 2 COMPLETE. Return only: ',return,'\nKeep all local caches; no further site pass is required.\n',sep='')
