# Coordinator consumes aggregate returns only. No source CLIF tables or site caches.
source('scripts/deployment_lib.R')
a<-commandArgs(trailingOnly=TRUE)
if(length(a)<2)stop('Usage: coordinator.R inspect|freeze CONFIG; coordinator.R predictors|collect CONFIG SITE_RETURN [...]')
clif_check_environment();command<-a[1];config<-deploy_read(a[2]);root<-normalizePath(config$output,mustWork=FALSE)
returns<-unlist(config$stage1_returns)
headers<-lapply(returns,function(p)deploy_read(file.path(p,'return_manifest.json')));sites<-vapply(headers,`[[`,character(1),'site')
if(length(sites)<3||anyDuplicated(sites))stop('At least three distinct sites are required for full and leave-one-site-out discovery.')
rosters<-setNames(lapply(DEPLOY_VARIANTS,function(v)unlist(config$discovery_sites[[v]]%||%as.list(sites))),DEPLOY_VARIANTS)
for(v in DEPLOY_VARIANTS)if(length(rosters[[v]])<3||anyDuplicated(rosters[[v]])||!setequal(rosters[[v]],sites))stop('Every participating site must contribute to every analysis. Repair the required discovery inputs; roster mismatch: ',v)
meta<-lapply(seq_along(returns),function(i){vs<-DEPLOY_VARIANTS[vapply(rosters,function(s)sites[i]%in%s,logical(1))];deploy_stage1_input(returns[i],config$accepted_stage1[[sites[i]]],vs)})
acceptance<-setNames(lapply(seq_along(sites),function(i)config$accepted_stage1[[sites[i]]]%||%list(site=sites[i],decision='same_source_contract',justification='Exact current source/environment contract.',return_manifest_md5=unname(tools::md5sum(file.path(returns[i],'return_manifest.json'))),source_hash=meta[[i]]$source_hash,r_lock_md5=meta[[i]]$r_lock_md5)),sites)
if(!(config$fingerprint_source%||%'exact_counts')%in%c('exact_counts','interval_reference'))stop('Unsupported fingerprint_source.')
Sys.setenv(CLIF_FINGERPRINT_SOURCE=config$fingerprint_source%||%'exact_counts',CLIF_OUTCOME_SITES=paste(sites,collapse=','),CLIF_ACCEPTED_LEGACY_INPUTS=if(any(vapply(acceptance,function(x)x$decision=='accepted_legacy_discovery',logical(1))))'true' else 'false')
params<-lapply(meta,function(m){p<-m$parameters;p$run_id<-NULL;p})
if(length(unique(vapply(params,jsonlite::toJSON,character(1),auto_unbox=TRUE)))!=1)stop('Sites used different scientific parameters.')
Sys.setenv(CLIF_FEDERATED_SITES=paste(sites,collapse=','),CLIF_APACHE_SITES=paste(sites,collapse=','),CLIF_REVISION_OUTPUT=file.path(root,'revision'),CLIF_REVISION_PRIMARY_INPUT=file.path(root,'primary'),CLIF_APACHE_OUTPUT=file.path(root,'apache'),CLIF_APACHE_PRIMARY=file.path(root,'primary'),R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
step<-function(key,script,args=character())deploy_step(root,key,script,args)
if(command=='inspect'){
 cat('Compatible Stage 1 returns from: ',paste(sites,collapse=', '),'\n')
 for(i in seq_along(returns)){
  z<-readr::read_csv(file.path(returns[i],'discovery/primary/silhouette.csv'),show_col_types=FALSE);best<-z[which.max(z$silhouette),]
  cat(sites[i],': local k=',best$k,'; silhouette=',round(best$silhouette,3),'\n',sep='')
 }
 cat('Review discovery/*/{clusters,fingerprints,silhouette}.csv, support missingness, cohort audits and provisional descriptors.\nEdit selected_k in coordinator config if desired; null uses the supported silhouette maximum. No outcomes are used to learn templates.\n');quit(status=0)
}
if(!command%in%c('freeze','predictors','collect'))stop('Unknown coordinator command.')
# Immutable coordination configuration prevents silently changing parameters on resume.
coordstate<-file.path(root,'local/coordinator_config.json')
if(file.exists(coordstate)){
 old<-deploy_read(coordstate);if(jsonlite::toJSON(old,auto_unbox=TRUE)!=jsonlite::toJSON(config,auto_unbox=TRUE))stop('Coordinator configuration changed. Choose a fresh output root.')
}else deploy_json(config,coordstate)
if(command=='freeze'){
 bundle<-file.path(root,'stage2_parameters')
 if(file.exists(file.path(bundle,'bundle_manifest.json'))){deploy_validate_inventory(bundle,'bundle_manifest.json');cat('Already frozen: ',bundle,'\n');quit(status=0)}
 for(v in DEPLOY_VARIANTS){
  work<-deploy_analysis_root(root,v);Sys.setenv(CLIF_FEDERATED_OUTPUT=work,CLIF_REVISION_VARIANT=v)
  selected<-config$selected_k[[v]]
  if(is.null(selected))Sys.unsetenv('CLIF_CONSENSUS_K') else Sys.setenv(CLIF_CONSENSUS_K=selected)
  Sys.setenv(CLIF_FEDERATED_SITES=paste(rosters[[v]],collapse=','))
  for(i in which(sites%in%rosters[[v]]))for(f in c('fingerprints.csv','clusters.csv','silhouette.csv','manifest.json'))deploy_copy(file.path(returns[i],'discovery',v,f),file.path(work,'exchange/discovery',sites[i],f))
  script<-if(v=='primary')'diagnostics/federated/central_discover.R' else if(v=='absorbing')'diagnostics/revision20261022/central_discover_absorbing.R' else 'diagnostics/revision20261022/central_discover.R'
  for(hold in c('none',sites))step(paste('freeze',v,hold,sep='_'),script,hold)
  deploy_tree(file.path(work,'coordinator'),file.path(bundle,v,'coordinator'))
 }
 primary<-deploy_read(file.path(root,'primary/coordinator/full/frozen_manifest.json'))
 deploy_json(list(stage='stage2_parameters',source_hash=deploy_hash(),r_lock_md5=deploy_lock_hash(),repo_commit=system2(Sys.which('git'),c('rev-parse','HEAD'),stdout=TRUE),sites=sites,variants=DEPLOY_VARIANTS,workspace_id=config$release_id%||%substr(primary$consensus_id,1,12),discovery_sites=rosters,stage1_acceptance=acceptance,fingerprint_source=config$fingerprint_source%||%'exact_counts',primary_consensus_id=primary$consensus_id,stage1_manifest_md5=setNames(as.list(vapply(returns,function(p)unname(tools::md5sum(file.path(p,'return_manifest.json'))),character(1))),sites),selected_k=config$selected_k,assignment='closest template for all support-sufficient stays; I retained for insufficient support',files=deploy_inventory(bundle)),file.path(bundle,'bundle_manifest.json'))
 cat('Frozen. Send this entire coordinator-created bundle to every participating site:\n',bundle,'\nSites run stage2-fit using it.\n');quit(status=0)
}
if(length(a)<3)stop('Supply all Stage 2 aggregate return directories after the coordinator config.')
fitreturns<-a[-c(1,2)];fm<-lapply(fitreturns,deploy_validate_inventory);ss<-vapply(fm,`[[`,character(1),'site')
if(anyDuplicated(ss)||!setequal(ss,sites))stop('Exactly one Stage 2 return per participating site is required.')
fitreturns<-fitreturns[match(sites,ss)];fm<-fm[match(sites,ss)]
for(m in fm)if(m$source_hash!=deploy_hash())stop('Stage 2 source contract mismatch.')
if(command=='predictors'){
 bundle<-file.path(root,'stage2_evaluation')
 if(file.exists(file.path(bundle,'bundle_manifest.json'))){deploy_validate_inventory(bundle,'bundle_manifest.json');cat('Already frozen: ',bundle,'\n');quit(status=0)}
 pm<-deploy_validate_inventory(file.path(root,'stage2_parameters'),'bundle_manifest.json')
 phash<-unname(tools::md5sum(file.path(root,'stage2_parameters/bundle_manifest.json')))
 for(m in fm)if(m$stage!='stage2_fit'||m$parameter_bundle_md5!=phash)stop('Return was not fit with the issued Stage 2 bundle.')
 for(v in setdiff(DEPLOY_VARIANTS,'absorbing')){
  work<-deploy_analysis_root(root,v);Sys.setenv(CLIF_FEDERATED_OUTPUT=work,CLIF_REVISION_VARIANT=v)
  for(i in seq_along(sites))deploy_tree(file.path(fitreturns[i],v),file.path(work,'sites',sites[i],'export'),c('csv','json'))
  step(paste0('stage_outcomes_',v),'diagnostics/federated/stage_outcomes.R')
  for(version in c('full',paste0('loo_',sites)))step(paste('pool',v,version,sep='_'),if(v=='primary')'diagnostics/federated/central_pool.R' else 'diagnostics/revision20261022/central_pool.R',version)
  if(v=='primary')step('pool_comparators','diagnostics/federated/central_comparators.R')
  deploy_tree(file.path(work,'coordinator'),file.path(bundle,v,'coordinator'))
 }
 # Absorbing groups are never fitted to the outcomes encoded in their inputs.
 deploy_tree(file.path(root,'stage2_parameters/absorbing/coordinator'),file.path(bundle,'absorbing/coordinator'))
 for(i in seq_along(sites)){
  deploy_tree(file.path(fitreturns[i],'review'),file.path(root,'revision/sites',sites[i],'export'),c('csv','json'))
  deploy_tree(file.path(fitreturns[i],'apache'),file.path(root,'apache/sites',sites[i],'export'),c('csv','json'))
 }
 Sys.setenv(CLIF_FEDERATED_OUTPUT=file.path(root,'primary'))
 step('freeze_scores','diagnostics/revision20261022/central_scores.R');step('freeze_apache','diagnostics/apache_adapted/central_train.R')
 step('shared_survival_grid','diagnostics/pooled_survival/central_grid.R',c(file.path(root,'survival_grid'),file.path(fitreturns,'survival_times')))
 deploy_tree(file.path(root,'revision/coordinator'),file.path(bundle,'scores/coordinator'))
 deploy_tree(file.path(root,'apache/coordinator'),file.path(bundle,'apache/coordinator'))
 deploy_tree(file.path(root,'survival_grid'),file.path(bundle,'survival_grid'))
 deploy_json(list(stage='stage2_evaluation',source_hash=deploy_hash(),r_lock_md5=pm$r_lock_md5,sites=sites,variants=pm$variants,workspace_id=pm$workspace_id,stage1_acceptance=pm$stage1_acceptance,primary_consensus_id=pm$primary_consensus_id,fit_manifest_md5=setNames(as.list(vapply(fitreturns,function(p)unname(tools::md5sum(file.path(p,'return_manifest.json'))),character(1))),sites),files=deploy_inventory(bundle)),file.path(bundle,'bundle_manifest.json'))
 cat('Send the frozen predictor/shared-grid bundle to all sites:\n',bundle,'\nSites run stage2-evaluate. No patient rows were read centrally.\n');quit(status=0)
}
for(m in fm)if(m$stage!='stage2_final')stop('collect requires final Stage 2 returns.')
Sys.setenv(CLIF_REQUIRE_FOUR_ENDPOINT_COX='true',CLIF_REQUIRE_POOLED_SURVIVAL='true')
step('final_collection','diagnostics/site_return/collect.R',c(file.path(root,'reports/site_returns'),file.path(fitreturns,'share')))
step('pooled_AJ_KM','diagnostics/pooled_survival/central_pool.R',c(file.path(root,'reports/pooled_survival'),file.path(root,'survival_grid'),file.path(fitreturns,'share/federated/pooled_survival/risksets')))
cat('Final site returns verified and competing-risk counts pooled in ',file.path(root,'reports'),'.\n')

python<-config$python %||% file.path(clif_root,'.venv',if(.Platform$OS.type=='windows')'Scripts/python.exe' else 'bin/python')
if(!file.exists(python))stop('Set coordinator config python to the plotting virtual environment, then rerun collect; completed steps are preserved.')
deploy_step(root,'pool_ROC','diagnostics/pooled_roc/central_plot.py',c(file.path(root,'reports/pooled_roc'),file.path(fitreturns,'pooled_roc')),python)
deploy_step(root,'plot_competing','diagnostics/pooled_survival/plot_competing.py',file.path(root,'reports/pooled_survival'),python)
