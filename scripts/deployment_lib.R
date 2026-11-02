source('scripts/runtime.R')
source('code/lib/load.R')
DEPLOY_VARIANTS<-c('primary','h24','h48','onset','support_only','low_flow','carry0','carry12','raw_support','absorbing')
deploy_hash<-function() {
 # Hash shared source only; ignored historical copies must not change a site contract.
 git<-Sys.which('git');if(!nzchar(git))stop('Git is required for source-contract checks; use a Git checkout of the coordinator-issued commit.')
 files<-system2(git,c('ls-files','--','code','diagnostics','scripts','renv.lock','config/waterfall_requirements.txt'),stdout=TRUE,stderr=TRUE)
 if(!is.null(attr(files,'status')))stop('Unable to read the tracked source inventory. Run from a Git checkout.')
 files<-sort(files[grepl('^(code|diagnostics|scripts)/.*[.](R|py|sh|bat)$',files)|files%in%c('renv.lock','config/waterfall_requirements.txt')])
 if(!length(files)||any(!file.exists(files)))stop('Tracked source files are missing; restore the coordinator-issued checkout.')
 digest::digest(as.list(setNames(unname(tools::md5sum(files)),files)),algo='sha256')
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
 if(!identical(cfg$analysis[names(expected)[names(expected)!='run_id']],expected[names(expected)[names(expected)!='run_id']]))stop('Scientific Stage 1 parameters differ from config/config_template.json. Only the coordinator may issue a new study specification.')
 cfg
}
deploy_env<-function(cfg,root,sites=cfg$site_name){
 env<-c(CLIF_RUN_CONFIG=normalizePath(Sys.getenv('CLIF_RUN_CONFIG')),CLIF_FEDERATED_OUTPUT=file.path(root,'primary'),CLIF_REVISION_OUTPUT=file.path(root,'revision'),CLIF_REVISION_PRIMARY_INPUT=file.path(root,'primary'),CLIF_APACHE_OUTPUT=file.path(root,'apache'),CLIF_APACHE_PRIMARY=file.path(root,'primary'),CLIF_WATERFALL_OUTPUT=file.path(root,'waterfall'),CLIF_FEDERATED_SITES=paste(sites,collapse=','),CLIF_APACHE_SITES=paste(sites,collapse=','),CLIF_LOCAL_DEPLOYMENT='true',R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep))
 env[paste0('CLIF_FEDERATED_INPUT_',cfg$site_name)]<-file.path(root,'waterfall')
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
 setNames(as.list(unname(tools::md5sum(files))),substring(files,nchar(dir)+2))
}
deploy_validate_inventory<-function(dir,manifest='return_manifest.json'){
 m<-deploy_read(file.path(dir,manifest));expected<-m$files;actual<-deploy_inventory(dir)
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
 deploy_json(c(list(status='complete',stage=stage,site=site,source_hash=deploy_hash(),r_version=as.character(getRversion()),r_lock_md5=unname(tools::md5sum('renv.lock')),files=deploy_inventory(dir)),extra),file.path(dir,'return_manifest.json'))
}
deploy_analysis_root<-function(root,variant)if(variant=='primary')file.path(root,'primary') else file.path(root,'revision/sensitivities',if(variant=='absorbing')'absorbing_v2' else variant)
