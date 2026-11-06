# Sourced only by scripts/site.R after environment, bundle and source validation.
if(!identical(bm$workflow,'one_pass_local_validation_v1'))stop('This bundle does not implement one-pass Stage2.')
for(v in DEPLOY_VARIANTS){
 dest<-file.path(deploy_analysis_root(root,v),'coordinator')
 deploy_tree(file.path(bundle,v,'coordinator'),dest,c('csv','json'))
 for(version in c('full',paste0('loo_',sites))){p<-deploy_load_prototypes(file.path(dest,version));saveRDS(p,file.path(dest,version,'prototypes.rds'))}
}
step('rebuild_raw','diagnostics/stage2/rebuild_inputs.R','raw')
step('rebuild_waterfall_bounds','diagnostics/waterfall/prepare.R')
step('rebuild_waterfall_test','diagnostics/waterfall/test_boundaries.py',python=cfg$python)
step('rebuild_waterfall','diagnostics/waterfall/run_waterfall.py',python=cfg$python)
step('rebuild_hourly','diagnostics/stage2/rebuild_inputs.R','waterfall')
for(v in DEPLOY_VARIANTS){
 Sys.setenv(CLIF_FEDERATED_OUTPUT=deploy_analysis_root(root,v),CLIF_REVISION_VARIANT=v)
 step(paste0('classify_',v),'diagnostics/stage2/classify_inputs.R',c(site,v))
}
step('finalize_rebuild','diagnostics/stage2/finalize_inputs.R',c(site,root,bundle))
Sys.setenv(CLIF_SUPPLEMENT_DERIVATION=file.path(cfg$work_root,'local/stage2_derivation.rds'))
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
invisible(deploy_stage2_context(cfg,bundle,bm,FALSE))
deploy_finish(return,'stage2',site,list(input_mode='rebuild_from_clif',stage1_outputs_required=FALSE,clustering_performed=FALSE,workflow=bm$workflow,primary_consensus_id=bm$primary_consensus_id,parameter_bundle_md5=unname(tools::md5sum(file.path(bundle,'bundle_manifest.json'))),validation='patient_grouped_local_cross_validation',external_site_outcome_validation=FALSE))
cat('STAGE 2 COMPLETE. Return only: ',return,'\nNo Stage1 results were required; no further site pass is required.\n',sep='')
