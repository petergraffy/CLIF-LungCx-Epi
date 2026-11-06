source('scripts/runtime.R')
source('code/lib/load.R')
DEPLOY_VARIANTS<-c('primary','h24','h48','onset','support_only','low_flow','carry0','carry12','raw_support','absorbing')
deploy_source_files<-function() {
 # Hash shared source only; ignored historical copies must not change a site contract.
 git<-Sys.which('git');if(!nzchar(git))stop('Git is required for source-contract checks; use a Git checkout of the coordinator-issued commit.')
 files<-system2(git,c('ls-files','--','code','diagnostics','scripts','renv.lock','config/waterfall_requirements.txt'),stdout=TRUE,stderr=TRUE)
 if(!is.null(attr(files,'status')))stop('Unable to read the tracked source inventory. Run from a Git checkout.')
 files<-sort(files[grepl('^(code|diagnostics|scripts)/.*[.](R|py|sh|bat)$',files)|files%in%c('renv.lock','config/waterfall_requirements.txt')],method='radix')
 if(!length(files)||any(!file.exists(files)))stop('Tracked source files are missing; restore the coordinator-issued checkout.')
 files
}
deploy_source_inventory<-function()as.list(source_text_md5(deploy_source_files()))
deploy_hash<-function()digest::digest(deploy_source_inventory(),algo='sha256')
deploy_lock_hash<-function()unname(source_text_md5('renv.lock'))
deploy_assert_contract<-function(m){
 if(!identical(m$source_hash,deploy_hash())){
  expected<-deploy_source_inventory();received<-m$source_files_md5
  # Older manifests serialized the same named hashes in locale-dependent order.
  # Accept that case only when its own hash reproduces AND every source file matches.
  ordered_match<-!is.null(received)&&setequal(names(received),names(expected))&&
   identical(received[names(expected)],expected)&&
   identical(m$source_hash,digest::digest(received,algo='sha256'))
  if(ordered_match){
   if(!identical(m$r_lock_md5,deploy_lock_hash()))stop('R lockfile contract mismatch.')
   return(invisible(TRUE))
  }
  changed<-if(is.null(received))'per-file provenance unavailable in this older return' else paste(sort(union(setdiff(names(received),names(expected)),names(expected)[vapply(names(expected),function(k)!identical(received[[k]],expected[[k]]),logical(1))])),collapse=', ')
  stop('Source contract mismatch: ',changed,'. Preserve the original return; use the coordinator-issued checkout and a fresh run_id. See code/README.md.')
 }
 if(!identical(m$r_lock_md5,deploy_lock_hash()))stop('R lockfile contract mismatch. Use the coordinator-issued lockfile; older byte-hashed returns need a coordinated refresh.')
 invisible(TRUE)
}
deploy_json<-function(x,path){dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);jsonlite::write_json(x,path,auto_unbox=TRUE,pretty=TRUE,null='null',digits=16)}
deploy_read<-function(path)jsonlite::fromJSON(path,simplifyVector=FALSE)
deploy_root<-function(cfg)file.path(cfg$repo,'output',paste0('run_',cfg$site_name,'_',cfg$analysis$run_id,'_federated'))
deploy_config<-function(path){
 cfg<-read_run_config(path);cfg$repo<-normalizePath(cfg$repo,winslash='/')
 if(cfg$repo!=clif_root)stop('config repo must match the current repository checkout.')
 if(cfg$file_type!='parquet')stop('The adopted waterfall deployment requires CLIF parquet tables.')
 if(is.null(cfg$python)||!file.exists(cfg$python))stop('Set python to the Python 3.12 virtual-environment executable in config/config.json.')
 expected<-read_run_config('config/config_template.json')$analysis
 keys<-setdiff(names(expected),'run_id')
 if(any(!vapply(keys,function(k)isTRUE(all.equal(cfg$analysis[[k]],expected[[k]],check.attributes=FALSE)),logical(1))))stop('Scientific analysis settings differ from config/config_template.json. Only the coordinator may issue a new study specification.')
 cfg
}
deploy_env<-function(cfg,root,sites=cfg$site_name){
 env<-c(CLIF_RUN_CONFIG=normalizePath(Sys.getenv('CLIF_RUN_CONFIG')),CLIF_FEDERATED_OUTPUT=file.path(root,'primary'),CLIF_REVISION_OUTPUT=file.path(root,'revision'),CLIF_REVISION_PRIMARY_INPUT=file.path(root,'primary'),CLIF_APACHE_OUTPUT=file.path(root,'apache'),CLIF_APACHE_PRIMARY=file.path(root,'primary'),CLIF_WATERFALL_OUTPUT=file.path(root,'waterfall'),CLIF_FEDERATED_SITES=paste(sites,collapse=','),CLIF_APACHE_SITES=paste(sites,collapse=','),CLIF_LOCAL_DEPLOYMENT='true',R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
 if(isTRUE(cfg$stage2_independent)){env['CLIF_STAGE2_NO_CLUSTERING']<-'true';env['CLIF_SOURCE_TIMEZONE']<-cfg$source_timezone}
 env[paste0('CLIF_FEDERATED_INPUT_',cfg$site_name)]<-file.path(cfg$stage1_root%||%root,'waterfall')
 env[paste0('CLIF_SITE_TABLES_',cfg$site_name)]<-cfg$tables_path
 env[paste0('CLIF_APACHE_TABLES_',cfg$site_name)]<-cfg$tables_path
 do.call(Sys.setenv,as.list(env));invisible(env)
}
deploy_step<-function(root,key,script,args=character(),python=NULL){
 dir<-file.path(root,'local/steps');dir.create(dir,recursive=TRUE,showWarnings=FALSE)
 marker<-file.path(dir,paste0(key,'.json'));log<-file.path(dir,paste0(key,'.log'))
 if(file.exists(marker)){message('Already complete: ',key);return(invisible(TRUE))}
 message('Running ',key,'; local log: ',log)
 command<-if(is.null(python))file.path(R.home('bin'),'Rscript') else python
 argv<-c(if(is.null(python))'--vanilla',shQuote(script),shQuote(args))
 status<-system2(command,argv,stdout=log,stderr=log)
 if(status!=0)stop(key,' failed. Inspect ',log,'. No complete return is declared. A partially written worker may require a fresh run_id.')
 deploy_json(list(status='passed',script=script),marker)
}
deploy_copy<-function(from,to){
 if(!file.exists(from))stop('Required input missing: ',from)
 dir.create(dirname(to),recursive=TRUE,showWarnings=FALSE)
 if(file.exists(to)){
  if(tools::md5sum(from)!=tools::md5sum(to))stop('Refusing to overwrite a different input: ',to)
 }else if(!file.copy(from,to))stop('Copy failed: ',to)
}
deploy_tree<-function(from,to,extensions=c('csv','json','rds')){
 files<-list.files(from,recursive=TRUE,full.names=TRUE)
 for(f in files){if(!tolower(tools::file_ext(f))%in%extensions)stop('Unexpected bundle file type: ',f);deploy_copy(f,file.path(to,substring(f,nchar(from)+2)))}
}
deploy_inventory<-function(dir) {
 files<-sort(list.files(dir,recursive=TRUE,full.names=TRUE))
 files<-files[!basename(files)%in%c('return_manifest.json','bundle_manifest.json')]
 if(any(file.info(files)$size==0))stop('Empty return/bundle file: ',paste(files[file.info(files)$size==0],collapse=', '))
 setNames(as.list(unname(tools::md5sum(files))),substring(files,nchar(dir)+2))
}
deploy_validate_inventory<-function(dir,manifest='return_manifest.json'){
 m<-deploy_read(file.path(dir,manifest));expected<-m$files;actual<-deploy_inventory(dir);actual[[manifest]]<-NULL
 if(!identical(sort(names(actual)),sort(names(expected))))stop('Return/bundle inventory differs from its manifest.')
 for(k in names(expected))if(!identical(actual[[k]],expected[[k]]))stop('Hash mismatch: ',k)
 m
}
deploy_check_csvs<-function(dir){
 forbidden<-c('patient_id','hospitalization_id','recorded_dttm','admission_dttm','discharge_dttm','t0','event_time','sequence')
 for(p in list.files(dir,pattern='[.]csv$',recursive=TRUE,full.names=TRUE)){
  cols<-names(readr::read_csv(p,n_max=0,show_col_types=FALSE));if(any(tolower(cols)%in%forbidden))stop('Patient-level field in return: ',p)
 }
 if(any(!tolower(tools::file_ext(list.files(dir,recursive=TRUE)))%in%c('csv','json','txt')))stop('Only aggregate CSV/JSON/TXT files may enter a site return.')
}
deploy_finish<-function(dir,stage,site,extra=list()){
 deploy_check_csvs(dir)
 deploy_json(c(list(status='complete',stage=stage,site=site,source_hash=deploy_hash(),source_hash_policy='text-md5-crlf-to-lf-v1',source_files_md5=deploy_source_inventory(),r_version=as.character(getRversion()),r_lock_md5=deploy_lock_hash(),files=deploy_inventory(dir)),extra),file.path(dir,'return_manifest.json'))
}
deploy_analysis_root<-function(root,variant)if(variant=='primary')file.path(root,'primary') else file.path(root,'revision/sensitivities',if(variant=='absorbing')'absorbing_v2' else variant)

# Explicitly scoped approvals bind original manifests; source hashes are never rewritten.
deploy_stage1_input<-function(dir,approval=NULL,variants=DEPLOY_VARIANTS){
 if(is.null(approval)){m<-deploy_validate_inventory(dir);if(m$stage!='stage1'||m$status!='complete')stop('Expected a completed Stage1 return.');deploy_assert_contract(m);return(m)}
 m<-deploy_read(file.path(dir,'return_manifest.json'))
 if(!approval$decision%in%c('accepted_legacy_discovery','same_source_contract')||!nzchar(approval$justification))stop('A documented discovery-snapshot acceptance is required.')
 if(m$stage!='stage1'||m$status!='complete'||m$site!=approval$site||m$source_hash!=approval$source_hash||m$r_lock_md5!=approval$r_lock_md5||
    unname(tools::md5sum(file.path(dir,'return_manifest.json')))!=approval$return_manifest_md5)stop('Original Stage1 snapshot differs from its approved lineage.')
 for(v in variants)for(name in c('fingerprints.csv','clusters.csv','silhouette.csv','manifest.json')){
  rel<-file.path('discovery',v,name);p<-file.path(dir,rel);expected<-m$files[[rel]]
  if(is.null(expected)||!file.exists(p)||file.info(p)$size==0||unname(tools::md5sum(p))!=expected)stop('Required discovery input failed integrity: ',m$site,'/',rel)
  if(grepl('[.]csv$',name)){cols<-names(readr::read_csv(p,n_max=0,show_col_types=FALSE));if(any(tolower(cols)%in%c('patient_id','hospitalization_id','t0','event_time','recorded_dttm','admission_dttm','discharge_dttm','sequence')))stop('Patient-level field in discovery input: ',rel)}
 }
 m
}
deploy_stage2_context<-function(cfg,bundle,bm,copy_caches=FALSE){
 site<-cfg$site_name
 if(bm$source_hash!=deploy_hash()||bm$r_lock_md5!=deploy_lock_hash())stop('Use the coordinator-issued Stage2 checkout and pinned environment.')
 if(!identical(bm$input_mode,'rebuild_from_clif'))stop('Use the independent Stage2 bundle shipped with this checkout.')
 if(!site%in%unlist(bm$sites))stop('Site not admitted to this release; contact the coordinator.')
 workspace_id<-bm$workspace_id
 if(!grepl('^[A-Za-z0-9_-]+$',workspace_id))stop('Invalid workspace ID.')
 root<-file.path(cfg$repo,'output',paste0('run_',site,'_',cfg$analysis$run_id,'_stage2_',workspace_id))
 variants<-unlist(bm$variants)
 for(v in variants)for(version in c('full',paste0('loo_',unlist(bm$sites)))){
  modeldir<-file.path(bundle,v,'coordinator',version);mp<-file.path(modeldir,'frozen_manifest.json');pp<-file.path(modeldir,'prototypes.rds')
  if(!file.exists(mp)||(!file.exists(pp)&&!file.exists(file.path(modeldir,'prototypes.csv'))))stop('Released templates are missing: ',v,'/',version)
  fm<-deploy_read(mp);proto<-deploy_load_prototypes(modeldir)
  if(!is.matrix(proto)||nrow(proto)!=fm$k||!ncol(proto)||any(!is.finite(proto))||any(proto<0)||anyDuplicated(rownames(proto))||anyDuplicated(colnames(proto)))stop('Invalid released prototype matrix: ',v,'/',version)
  if(v=='primary'&&nrow(proto)!=4L)stop('This release requires four primary phenotypes.')
  if(version!='full'&&sub('^loo_','',version)%in%unlist(fm$training_sites))stop('Held-out site leaked into template learning.')
  if(v=='primary'&&version=='full'&&fm$consensus_id!=bm$primary_consensus_id)stop('Primary consensus identity mismatch.')
 }

 cfg$stage1_root<-NULL;cfg$work_root<-file.path(root,'inputs');cfg$stage2_independent<-TRUE
 cfg$source_timezone<-bm$source_timestamp_policy[[site]]$source_timezone%||%'UTC'
 if(!cfg$source_timezone%in%c('UTC',OlsonNames()))stop('Invalid source timezone in frozen bundle.')
 files<-sort(list.files(cfg$tables_path,pattern='[.]parquet$',full.names=TRUE,recursive=TRUE,ignore.case=TRUE))
 info<-file.info(files)
 parameter_names<-sort(union(names(cfg$analysis),names(run_defaults())),method='radix')
 parameters<-setNames(lapply(parameter_names,function(k)cfg$analysis[[k]]),parameter_names)
 state<-list(site=site,source_hash=deploy_hash(),bundle_md5=unname(tools::md5sum(file.path(bundle,'bundle_manifest.json'))),parameters=parameters,tables_path=cfg$tables_path,source_timezone=cfg$source_timezone,files=setNames(lapply(seq_along(files),function(i)list(bytes=info$size[i],mtime=as.numeric(info$mtime[i]))),files))
 statefile<-file.path(root,'local/rebuild_inputs.json')
 if(file.exists(statefile)&&!identical(deploy_read(statefile),jsonlite::fromJSON(jsonlite::toJSON(state,auto_unbox=TRUE,digits=16,null='null'),simplifyVector=FALSE)))stop('Source tables, code, config or bundle changed since Stage2 started; use a fresh analysis.run_id.')
 config<-file.path(root,'local/runtime_config.json')
 if(copy_caches){
  deploy_json(state,statefile);deploy_json(cfg,config)
 }
 list(cfg=cfg,root=root,config=config,input_mode='rebuild_from_clif')
}

# Released templates are portable aggregate CSVs; RDS caches are created only locally.
deploy_load_prototypes<-function(dir){
 csv<-file.path(dir,'prototypes.csv')
 if(!file.exists(csv))return(readRDS(file.path(dir,'prototypes.rds')))
 z<-readr::read_csv(csv,show_col_types=FALSE);ph<-sort(unique(z$phenotype),method='radix');f<-unique(z$feature)
 if(anyDuplicated(z[c('phenotype','feature')])||nrow(z)!=length(ph)*length(f))stop('Incomplete prototype CSV.')
 p<-matrix(NA_real_,length(ph),length(f),dimnames=list(ph,f));p[cbind(match(z$phenotype,ph),match(z$feature,f))]<-z$probability;p
}

# Publish aggregate CSV/JSON templates only, without private coordinator paths or RDS caches.
deploy_release_templates<-function(from,to){
 for(f in list.files(from,recursive=TRUE,full.names=TRUE)){
  rel<-substring(f,nchar(from)+2L);if(!tools::file_ext(f)%in%c('csv','json'))next
  dest<-file.path(to,rel)
  if(tools::file_ext(f)=='json'){
   m<-deploy_read(f);m$discovery_inputs<-NULL
   if(!is.null(m$source_hashes))m$source_hashes<-lapply(m$source_hashes,function(x)setNames(x,basename(names(x))))
   deploy_json(m,dest)
  }else deploy_copy(f,dest)
 }
}
