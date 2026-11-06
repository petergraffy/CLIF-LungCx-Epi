source('scripts/deployment_lib.R');source('diagnostics/federated/contract.R')
base<-tempfile('run_stage2_workspace_',tmpdir='output');dir.create(base);base<-normalizePath(base)
site<-'SYNTHETIC_SITE';cfg<-list(repo=base,site_name=site,analysis=run_defaults());cfg$analysis$run_id<-'original'
original<-deploy_root(cfg);ret<-file.path(original,'stage1_return');dir.create(ret,recursive=TRUE)
t0<-utc('2020-01-01');ids<-c('synthetic_a','synthetic_b')
for(v in DEPLOY_VARIANTS){
 w<-if(v=='h24')24 else if(v=='h48')48 else 72
 data<-tibble::tibble(hospitalization_id=ids,t0=t0,time_from_landmark=10)
 panel<-tidyr::expand_grid(hospitalization_id=ids,h=0:(w-1))
 cache<-file.path(deploy_analysis_root(original,v),'sites',site,'local/discovery.rds');dir.create(dirname(cache),recursive=TRUE)
 saveRDS(list(data=data,panel=panel,input_root=base),cache)
 for(name in c('fingerprints.csv','clusters.csv','silhouette.csv','manifest.json')){
  rel<-file.path('discovery',v,name);dest<-file.path(ret,rel);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
  if(grepl('csv$',name))writeLines('synthetic,n\nvalue,2',dest) else writeLines('{"synthetic": true}',dest)
  export<-file.path(deploy_analysis_root(original,v),'sites',site,'export/discovery',name);dir.create(dirname(export),recursive=TRUE,showWarnings=FALSE);file.copy(dest,export)
 }
}
deploy_finish(ret,'stage1',site,list(parameters=cfg$analysis))
m<-deploy_read(file.path(ret,'return_manifest.json'));oldsource<-m$source_hash
# Represent an explicitly accepted original checkout; its hash is preserved unchanged.
m$source_hash<-'legacy_synthetic_source';deploy_json(m,file.path(ret,'return_manifest.json'))
deploy_json(list(source_hash=m$source_hash),file.path(original,'local/deployment.json'))
approval<-list(site=site,decision='accepted_legacy_discovery',justification='Synthetic migration test.',source_hash=m$source_hash,r_lock_md5=m$r_lock_md5,return_manifest_md5=unname(tools::md5sum(file.path(ret,'return_manifest.json'))))
bm<-list(source_hash=deploy_hash(),r_lock_md5=deploy_lock_hash(),sites=site,primary_consensus_id='synthetic_consensus',workspace_id='release_test',stage1_acceptance=setNames(list(approval),site),variants=DEPLOY_VARIANTS)
for(v in DEPLOY_VARIANTS)for(version in c('full',paste0('loo_',site))){
 path<-file.path(base,v,'coordinator',version);dir.create(path,recursive=TRUE)
 proto<-diag(4);dimnames(proto)<-list(paste0('P',1:4),paste0('feature',1:4));saveRDS(proto,file.path(path,'prototypes.rds'))
 deploy_json(list(k=4,training_sites=if(version=='full')list(site) else list(),consensus_id='synthetic_consensus'),file.path(path,'frozen_manifest.json'))
}
original_hash<-unname(tools::md5sum(file.path(ret,'return_manifest.json')))
x<-deploy_stage2_context(cfg,base,bm,FALSE);stopifnot(!dir.exists(x$root))
y<-deploy_stage2_context(cfg,base,bm,TRUE);stopifnot(y$root!=original,dir.exists(y$root))
stopifnot(unname(tools::md5sum(file.path(ret,'return_manifest.json')))==original_hash,deploy_read(file.path(original,'local/deployment.json'))$source_hash=='legacy_synthetic_source')
for(v in DEPLOY_VARIANTS){a<-file.path(deploy_analysis_root(original,v),'sites',site,'local/discovery.rds');b<-file.path(deploy_analysis_root(y$root,v),'sites',site,'local/discovery.rds');stopifnot(unname(tools::md5sum(a))==unname(tools::md5sum(b)))}
# Resuming is deterministic; copying caches never edits the original snapshot.
invisible(deploy_stage2_context(cfg,base,bm,TRUE))
bad<-bm;bad$sites<-'OTHER';stopifnot(inherits(tryCatch(deploy_stage2_context(cfg,base,bad,FALSE),error=identity),'error'))
# A fresh Stage2 checkout can reference caches in the preserved old checkout.
newcfg<-cfg;newcfg$repo<-file.path(base,'new_checkout');dir.create(newcfg$repo);newcfg$stage1_root<-original
z<-deploy_stage2_context(newcfg,base,bm,FALSE);stopifnot(z$original==original,startsWith(z$root,newcfg$repo))
cat('PASS: read-only Stage2 preflight, isolated cache migration, original hashes preserved, reproducible resume, explicit admission and separate checkout\n')
