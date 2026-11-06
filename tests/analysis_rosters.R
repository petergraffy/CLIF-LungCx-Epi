source('scripts/deployment_lib.R')
base<-tempfile('run_analysis_rosters_',tmpdir='output');dir.create(base);base<-normalizePath(base)
ids<-c('ALPHA','BETA_2','GAMMA-3','DELTA');Sys.setenv(CLIF_SYNTHETIC_SITES=paste(ids,collapse=','))
for(v in DEPLOY_VARIANTS){
 rc<-system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','tests/deployment_fingerprints.R',v,shQuote(base)),stdout=file.path(base,paste0(v,'.log')),stderr=file.path(base,paste0(v,'.log')));stopifnot(rc==0)
}
returns<-file.path(base,ids)
for(p in returns)deploy_finish(p,'stage1',basename(p),list(parameters=run_defaults()))
p<-returns[4];m<-deploy_read(file.path(p,'return_manifest.json'));m$source_hash<-'accepted_legacy_test'
# Ancillary descriptors do not enter discovery, but every site's discovery does.
m$files[['descriptives/pending.csv']]<-'pending';deploy_json(m,file.path(p,'return_manifest.json'))
approval<-list(site='DELTA',decision='accepted_legacy_discovery',justification='Synthetic scoped acceptance.',source_hash=m$source_hash,r_lock_md5=m$r_lock_md5,return_manifest_md5=unname(tools::md5sum(file.path(p,'return_manifest.json'))))
cfg<-list(output=file.path(base,'central'),stage1_returns=returns,selected_k=list(primary=4),accepted_stage1=list(DELTA=approval))
config<-file.path(base,'coordinator.json')
run<-function(label){deploy_json(cfg,config);system2(file.path(R.home('bin'),'Rscript'),c('--vanilla','scripts/coordinator.R','freeze',shQuote(config)),stdout=file.path(base,paste0(label,'.log')),stderr=file.path(base,paste0(label,'.log')))}
# Missing inputs are repair items. They must not silently remove a site.
h24<-file.path(p,'discovery/h24');backup<-file.path(base,'h24_backup');dir.create(backup);file.copy(list.files(h24,full.names=TRUE),backup);unlink(h24,recursive=TRUE)
stopifnot(run('missing')!=0)
dir.create(h24);file.copy(list.files(backup,full.names=TRUE),h24)
cfg$discovery_sites<-list(h24=as.list(ids[1:3]));stopifnot(run('partial_roster')!=0);cfg$discovery_sites<-NULL
rc<-run('freeze')
if(rc!=0){cat(tail(readLines(file.path(base,'freeze.log')),20),sep='\n');stop('All-site freeze failed')}
bundle<-file.path(base,'central/stage2_parameters');b<-deploy_validate_inventory(bundle,'bundle_manifest.json')
full<-deploy_read(file.path(bundle,'h24/coordinator/full/frozen_manifest.json'));held<-deploy_read(file.path(bundle,'h24/coordinator/loo_DELTA/frozen_manifest.json'))
stopifnot(setequal(unlist(full$training_sites),ids),setequal(unlist(full$outcome_training_sites),ids),!('DELTA'%in%unlist(held$training_sites)),!('DELTA'%in%unlist(held$outcome_training_sites)))
primary<-deploy_read(file.path(bundle,'primary/coordinator/full/frozen_manifest.json'));stopifnot(setequal(unlist(primary$training_sites),ids),primary$k==4,primary$fingerprint_source=='exact_counts')
write('tampered',file=file.path(p,'discovery/primary/fingerprints.csv'),append=TRUE)
stopifnot(inherits(tryCatch(deploy_stage1_input(p,approval,'primary'),error=identity),'error'))
cat('PASS: all sites contribute to every analysis; missing discovery files and partial rosters rejected; ancillary descriptors do not alter training; held-out validation and exact-count four-phenotype primary preserved\n')
