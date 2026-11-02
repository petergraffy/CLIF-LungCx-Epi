source('scripts/deployment_lib.R')
clif_check_environment()
check<-function(x,label){stopifnot(isTRUE(x));cat('PASS:',label,'\n')}
# No source CLIF data is accessed in this coordinator/contract round trip.
base<-tempfile('run_deployment_',tmpdir='output');dir.create(base,recursive=TRUE)
base<-normalizePath(base)
returns<-file.path(base,c('ALPHA','BETA_2','GAMMA-3'))
for(v in DEPLOY_VARIANTS){
 status<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','tests/deployment_fingerprints.R',v,shQuote(base)),stdout=file.path(base,paste0(v,'.log')),stderr=file.path(base,paste0(v,'.log')))
 check(status==0,paste('synthetic discovery',v))
}
for(i in seq_along(returns))deploy_finish(returns[i],'stage1',basename(returns[i]),list(variants=DEPLOY_VARIANTS,parameters=run_defaults(),shared_phenotypes_assigned=FALSE,local_outcome_models_fit=FALSE))
coordinator<-file.path(base,'coordinator.json');config<-list(output=file.path(base,'central'),stage1_returns=returns,selected_k=list(primary=4))
deploy_json(config,coordinator)
run<-function(action){log<-file.path(base,paste0(action,'.log'));status<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','scripts/coordinator.R',action,shQuote(coordinator)),stdout=log,stderr=log);if(status!=0)cat(tail(readLines(log),10),sep='\n');status}
check(run('inspect')==0,'coordinator inspect accepts arbitrary site IDs')
check(run('freeze')==0,'coordinator freezes full and held-out templates for all variants')
bundle<-file.path(base,'central/stage2_parameters');m<-deploy_validate_inventory(bundle,'bundle_manifest.json')
check(setequal(unlist(m$sites),basename(returns)),'bundle contains the actual participating sites')
for(v in DEPLOY_VARIANTS)for(hold in basename(returns)){
 p<-deploy_read(file.path(bundle,v,'coordinator',paste0('loo_',hold),'frozen_manifest.json'))
 check(!hold%in%unlist(p$training_sites),paste(v,hold,'excluded from held-out learning'))
}
p<-deploy_read(file.path(bundle,'primary/coordinator/full/frozen_manifest.json'))
check(p$k==4&&p$selection=='coordinator_selected_supported_k','coordinator parameter choice is recorded and respected')
check(run('freeze')==0,'unchanged completed freeze resumes without overwriting')
f<-file.path(returns[1],'discovery/primary/silhouette.csv');write('tampered',file=f,append=TRUE)
check(inherits(tryCatch(deploy_validate_inventory(returns[1]),error=identity),'error'),'tampered site return is rejected')
# The return gate rejects identifier-bearing CSVs, even with a valid directory shape.
bad<-file.path(base,'bad');dir.create(bad);readr::write_csv(tibble::tibble(hospitalization_id='synthetic'),file.path(bad,'forbidden.csv'))
check(inherits(tryCatch(deploy_finish(bad,'stage1','BAD'),error=identity),'error'),'patient-level CSV headers cannot be packaged')
check(!any(grepl('(^local/|[.]rds$)',names(m$files)))||all(grepl('prototypes[.]rds$',names(m$files)[grepl('[.]rds$',names(m$files))])),'only aggregate prototype RDS files are in coordinator bundle')
cat('Deployment synthetic round trip completed in ',base,'\n',sep='')
