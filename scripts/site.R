# Site-facing entry point. Run with Rscript --vanilla scripts/site.R ...
source('scripts/deployment_lib.R')
a<-commandArgs(trailingOnly=TRUE)
if(length(a)<2)stop('Usage: scripts/site.R preflight|stage1|stage2-preflight|stage2 CONFIG [COORDINATOR_BUNDLE]')
command<-a[1];config<-normalizePath(a[2],mustWork=TRUE);Sys.setenv(CLIF_RUN_CONFIG=config)
clif_check_environment();cfg<-deploy_config(config);site<-cfg$site_name;root<-deploy_root(cfg)
if(!command%in%c('preflight','stage1','stage2-preflight','stage2'))stop('Unknown command.')
# Source availability checks only: table contents are read by the local analysis.
files<-list.files(cfg$tables_path,pattern='[.]parquet$',recursive=TRUE,full.names=TRUE,ignore.case=TRUE)
bases<-tolower(sub('^clif_','',tools::file_path_sans_ext(basename(files)),ignore.case=TRUE))
required<-c('patient','hospitalization','adt','hospital_diagnosis','respiratory_support','vitals','labs','medication_admin_continuous','patient_assessments')
for(t in required)if(sum(bases==t)!=1)stop('Expected exactly one parquet table for ',t,'; found ',sum(bases==t))
status<-system2(cfg$python,c(shQuote('scripts/python_preflight.py')),stdout=TRUE,stderr=TRUE)
if(!is.null(attr(status,'status')))stop('Python preflight failed: ',paste(status,collapse='\n'))
if(command=='preflight'){
 cat('PASS: pinned R packages, Python waterfall and required source tables.\n')
 cat('Site: ',site,'\nStage 1 return: ',file.path(root,'stage1_return'), '\n',sep='')
 cat('This check does not derive cohorts or contact the coordinator.\n');quit(status=0)
}
if(command%in%c('stage2-preflight','stage2')){
 if(length(a)>3L)stop('Supply at most one bundle directory.')
 bundle<-normalizePath(if(length(a)==3L)a[3] else 'frozen/stage2_v1',mustWork=TRUE);bm<-deploy_validate_inventory(bundle,'bundle_manifest.json')
 if(command=='stage2-preflight'&&bm$stage!='stage2_parameters')stop('Preflight requires the Stage2 parameter bundle.')
 if(command=='stage2'&&bm$stage!='stage2_parameters')stop('Stage2 requires the released parameter bundle.')
 context<-deploy_stage2_context(cfg,bundle,bm,copy_caches=command!='stage2-preflight');cfg<-context$cfg;root<-context$root
 if(command=='stage2-preflight'){
  cat('PASS: Stage2 code/environment, approved Stage1 lineage, and all local discovery caches.\n')
  cat('Original Stage1 preserved: ',context$original,'\nStage2 workspace: ',root,'\n',sep='');quit(status=0)
 }
 dir.create(file.path(root,'local'),recursive=TRUE,showWarnings=FALSE)
}else{
dir.create(file.path(root,'local'),recursive=TRUE,showWarnings=FALSE)
statefile<-file.path(root,'local/deployment.json')
state<-list(site=site,parameters=cfg$analysis,source_hash=deploy_hash(),r_lock_md5=deploy_lock_hash())
if(file.exists(statefile)){
 old<-deploy_read(statefile)
 if(old$source_hash!=state$source_hash||!identical(jsonlite::toJSON(old$parameters,auto_unbox=TRUE),jsonlite::toJSON(state$parameters,auto_unbox=TRUE)))stop('Code/environment/settings changed since Stage 1. Use the original checkout or a new coordinated run_id.')
}else{
 if(command!='stage1')stop('Run Stage 1 locally before Stage 2.')
 deploy_json(state,statefile)
}
}
deploy_env(cfg,root)
step<-function(key,script,args=character(),python=NULL)deploy_step(root,key,script,args,python)
if(command=='stage1'){
 return<-file.path(root,'stage1_return')
 if(file.exists(file.path(return,'return_manifest.json'))){deploy_validate_inventory(return);cat('Stage 1 already complete. Return only: ',return,'\n');quit(status=0)}
 step('01_prepare','scripts/prepare_site.R')
 step('02_waterfall_bounds','diagnostics/waterfall/prepare.R')
 step('03_waterfall_test','diagnostics/waterfall/test_boundaries.py',python=cfg$python)
 step('04_waterfall','diagnostics/waterfall/run_waterfall.py',python=cfg$python)
 step('05_waterfall_derive','diagnostics/waterfall/compare.R')
 step('06_waterfall_QA','diagnostics/waterfall/qa.R')
 for(v in DEPLOY_VARIANTS){
  Sys.setenv(CLIF_FEDERATED_OUTPUT=deploy_analysis_root(root,v),CLIF_REVISION_VARIANT=v)
  script<-if(v=='primary')'diagnostics/federated/site_discover.R' else if(v=='absorbing')'diagnostics/revision20261022/site_discover_absorbing.R' else 'diagnostics/revision20261022/site_discover.R'
  step(paste0('discovery_',v),script,site)
  from<-file.path(deploy_analysis_root(root,v),'sites',site,'export/discovery')
  for(f in c('fingerprints.csv','clusters.csv','silhouette.csv','manifest.json'))deploy_copy(file.path(from,f),file.path(return,'discovery',v,f))
 }
 Sys.setenv(CLIF_FEDERATED_OUTPUT=file.path(root,'primary'))
 step('07_descriptors_supplement','scripts/stage1_descriptors.R')
 deploy_tree(file.path(root,'stage1_descriptors/share/descriptives'),file.path(return,'descriptives'),c('csv','json'))
 baseline<-init_run(cfg)
 deploy_tree(file.path(baseline,'share/supplement'),file.path(return,'supplement'),c('csv','json','txt'))
 for(f in c('cohort_exclusions.csv','cohort_by_year.csv','cohort_patient_counts.csv','outcome_timestamp_audit.csv','source_availability.csv','analysis_plan.csv','arf_detection_audit.csv','derivation_status.csv'))deploy_copy(file.path(baseline,'share',f),file.path(return,'cohort_audit',f))
 writeLines(c('STAGE 1 COMPLETE. Return this folder only.','No shared phenotype has been selected or assigned.','Discovery inputs: discovery/. Descriptors and supplement are review diagnostics, not consensus features.','Preserve this checkout, the site config and all output/run_SITE_RUN_ID* folders locally.','STOP and wait for the coordinator Stage 2 bundle. Do not run legacy consensus/model commands.'),file.path(return,'READ_ME_FIRST.txt'))
 deploy_copy('renv.lock',file.path(return,'environment/r_lock.json'))
 deploy_copy('config/waterfall_requirements.txt',file.path(return,'environment/waterfall_requirements.txt'))
 writeLines(capture.output(sessionInfo()),file.path(return,'environment/R_sessionInfo.txt'))
 deploy_finish(return,'stage1',site,list(variants=DEPLOY_VARIANTS,parameters=cfg$analysis,shared_phenotypes_assigned=FALSE,local_outcome_models_fit=FALSE,python_environment=status))
 cat('\nSTAGE 1 COMPLETE. STOP HERE.\nReturn only: ',return,'\nKeep all local caches. Wait for the coordinator Stage 2 bundle.\n',sep='');quit(status=0)
}
sites<-unlist(bm$sites);if(!site%in%sites)stop('This site is not included in the coordinator bundle.')
deploy_env(cfg,root,sites)
source('scripts/stage2_run.R')
