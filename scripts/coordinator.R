# Coordinator consumes aggregate returns only. No source CLIF tables or site caches.
source('scripts/deployment_lib.R')
a<-commandArgs(trailingOnly=TRUE)
if(length(a)<2)stop('Usage: coordinator.R inspect|freeze CONFIG; coordinator.R collect CONFIG SITE_RETURN [...]')
clif_check_environment();command<-a[1];config<-deploy_read(a[2]);root<-normalizePath(config$output,mustWork=FALSE)
if(command=='collect'){source('scripts/stage2_collect.R');quit(status=0)}
returns<-unlist(config$stage1_returns)
headers<-lapply(returns,function(p)deploy_read(file.path(p,'return_manifest.json')));sites<-vapply(headers,`[[`,character(1),'site')
if(length(sites)<3||anyDuplicated(sites))stop('At least three distinct sites are required for full and leave-one-site-out discovery.')
rosters<-setNames(lapply(DEPLOY_VARIANTS,function(v)unlist(config$discovery_sites[[v]]%||%as.list(sites))),DEPLOY_VARIANTS)
for(v in DEPLOY_VARIANTS){
 roster<-rosters[[v]];omitted<-setdiff(sites,roster);exceptions<-config$discovery_exceptions[[v]]
 if(length(roster)<3||anyDuplicated(roster)||!all(roster%in%sites)||any(!vapply(omitted,function(site)!is.null(exceptions[[site]])&&nzchar(exceptions[[site]]),logical(1))))stop('Every discovery omission needs an explicit coordinator-approved exception: ',v)
}
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
if(!command%in%c('freeze','collect'))stop('Unknown coordinator command.')
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
  deploy_release_templates(file.path(work,'coordinator'),file.path(bundle,v,'coordinator'))
 }
 primary<-deploy_read(file.path(root,'primary/coordinator/full/frozen_manifest.json'))
 deploy_copy('config/stage2_analysis.json',file.path(bundle,'analysis_plan.json'))
 writeLines(c('# Stage 2: one run, one return', '', 'Follow [the complete site instructions](../../code/README.md).', '', 'Restore the pinned R/Python environments, configure the CLIF table paths; Stage 1 results are not needed, then run from the repository root:', '', '```sh', 'Rscript --vanilla scripts/site.R stage2-preflight config/config.json', 'Rscript --vanilla scripts/site.R stage2 config/config.json', '```', '', 'Return only the printed stage2_return/ folder. No second site pass is required.', 'The bundle contains aggregate templates and a fixed analysis plan; patient data stay local.'),file.path(bundle,'SITE_INSTRUCTIONS.md'))
 deploy_json(list(stage='stage2_parameters',workflow='one_pass_local_validation_v1',input_mode='rebuild_from_clif',source_timestamp_policy=config$source_timestamp_policy,stage1_reference_counts=config$stage1_reference_counts,discovery_exceptions=config$discovery_exceptions,source_hash=deploy_hash(),r_lock_md5=deploy_lock_hash(),repo_commit=system2(Sys.which('git'),c('rev-parse','HEAD'),stdout=TRUE),sites=sites,variants=DEPLOY_VARIANTS,workspace_id=config$release_id%||%substr(primary$consensus_id,1,12),discovery_sites=rosters,stage1_acceptance=acceptance,fingerprint_source=config$fingerprint_source%||%'exact_counts',primary_consensus_id=primary$consensus_id,stage1_manifest_md5=setNames(as.list(vapply(returns,function(p)unname(tools::md5sum(file.path(p,'return_manifest.json'))),character(1))),sites),selected_k=config$selected_k,assignment='closest template for all support-sufficient stays; I retained for insufficient support',files=deploy_inventory(bundle)),file.path(bundle,'bundle_manifest.json'))
 cat('Frozen. Send this entire coordinator-created bundle to every participating site:\n',bundle,'\nSites run scripts/site.R stage2 once using it.\n');quit(status=0)
}
