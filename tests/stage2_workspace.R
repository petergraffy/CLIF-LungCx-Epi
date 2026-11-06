source('scripts/deployment_lib.R')
base<-tempfile('run_stage2_workspace_',tmpdir='output');dir.create(base);base<-normalizePath(base)
site<-'SYNTHETIC_SITE';tables<-file.path(base,'tables');dir.create(tables);writeLines('fixture metadata only',file.path(tables,'synthetic.parquet'))
cfg<-list(repo=base,tables_path=tables,site_name=site,analysis=run_defaults(),stage1_root=file.path(base,'DELETED_STAGE1'));cfg$analysis$run_id<-'new_stage2_run'
bm<-list(source_hash=deploy_hash(),r_lock_md5=deploy_lock_hash(),sites=list(site),primary_consensus_id='synthetic_consensus',workspace_id='independent_test',input_mode='rebuild_from_clif',variants=as.list(DEPLOY_VARIANTS))
for(v in DEPLOY_VARIANTS)for(version in c('full',paste0('loo_',site))){
 p<-file.path(base,v,'coordinator',version);dir.create(p,recursive=TRUE)
 proto<-diag(4);dimnames(proto)<-list(paste0('P',1:4),paste0('feature',1:4));saveRDS(proto,file.path(p,'prototypes.rds'))
 deploy_json(list(k=4,training_sites=if(version=='full')list(site) else list(),consensus_id='synthetic_consensus'),file.path(p,'frozen_manifest.json'))
}
x<-deploy_stage2_context(cfg,base,bm,FALSE);stopifnot(!dir.exists(x$root),is.null(x$cfg$stage1_root),x$cfg$work_root==file.path(x$root,'inputs'))
y<-deploy_stage2_context(cfg,base,bm,TRUE);stopifnot(file.exists(y$config),!dir.exists(deploy_root(cfg)),!dir.exists(cfg$stage1_root))
invisible(deploy_stage2_context(cfg,base,bm,TRUE))
invisible(deploy_stage2_context(read_run_config(y$config),base,bm,FALSE))
bad<-bm;bad$sites<-list('OTHER');stopifnot(inherits(tryCatch(deploy_stage2_context(cfg,base,bad,FALSE),error=identity),'error'))
writeLines('changed source',file.path(tables,'synthetic.parquet'))
stopifnot(inherits(tryCatch(deploy_stage2_context(cfg,base,bm,FALSE),error=identity),'error'))
cfg$analysis$run_id<-'another_run';stopifnot(!dir.exists(deploy_stage2_context(cfg,base,bm,FALSE)$root))
# Source policies preserve aware instants and handle DST ambiguities as missing.
Sys.setenv(CLIF_SOURCE_TIMEZONE='America/Chicago')
x<-utc('2020-01-01 12:00:00');stopifnot(source_time(x,arrow::timestamp('us'))==utc('2020-01-01 18:00:00'),source_time(x,arrow::timestamp('us',timezone='UTC'))==x,is.na(source_time(utc('2020-11-01 01:30:00'),arrow::timestamp('us'))))
Sys.setenv(CLIF_STAGE2_NO_CLUSTERING='true')
stopifnot(inherits(tryCatch(sequence_fit(NULL),error=identity),'error'),inherits(tryCatch(bootstrap_stability(NULL,NULL,1,1),error=identity),'error'))
cat('PASS: Stage2 needs no Stage1 files; read-only preflight, isolated workspace, safe resume/source-change rejection, new run IDs, explicit admission, timestamp/DST policies and no-clustering guards\n')
