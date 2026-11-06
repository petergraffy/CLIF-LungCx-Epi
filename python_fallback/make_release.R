# Coordinator only: bind the I/O-only overrides to the unchanged frozen release.
source('scripts/deployment_lib.R')
parent<-deploy_validate_inventory('frozen/stage2_v1','bundle_manifest.json')
stopifnot(identical(deploy_hash(),parent$source_hash),identical(deploy_lock_hash(),parent$r_lock_md5))
original<-deploy_source_files();scratch<-tempfile('run_python_release_',tmpdir='output')
dir.create(scratch);scratch<-normalizePath(scratch)
for(p in original){dir.create(dirname(file.path(scratch,p)),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(p,file.path(scratch,p)))}
overrides<-list.files('python_fallback/overrides',full.names=TRUE,recursive=TRUE)
for(p in overrides){rel<-substring(p,nchar('python_fallback/overrides/')+1L);stopifnot(rel%in%original);stopifnot(file.copy(p,file.path(scratch,rel),overwrite=TRUE))}
dir.create(file.path(scratch,'fallback'))
for(n in c('parquet.R','parquet_bridge.py'))stopifnot(file.copy(file.path('python_fallback',n),file.path(scratch,'fallback',n)))
effective<-sort(c(original,'fallback/parquet.R','fallback/parquet_bridge.py'),method='radix')
old<-getwd();setwd(scratch);hashes<-as.list(source_text_md5(effective));setwd(old)
distribution<-sort(list.files('python_fallback',pattern='[.](R|py)$',recursive=TRUE,full.names=TRUE),method='radix')
release<-list(release='stage2_python_v1',parent_bundle_md5=unname(tools::md5sum('frozen/stage2_v1/bundle_manifest.json')),parent_source_hash=parent$source_hash,r_lock_md5=parent$r_lock_md5,parent_source_files_md5=deploy_source_inventory(),parent_aux_files_md5=as.list(source_text_md5(c('config/config_template.json','config/stage2_analysis.json'))),effective_source_hash=digest::digest(hashes,algo='sha256'),effective_source_files_md5=hashes,distribution_files_md5=as.list(source_text_md5(distribution)),scientific_changes=FALSE,r_arrow_required=FALSE)
deploy_json(release,'python_fallback/release.json')
python<-Sys.getenv('CLIF_TEST_WATERFALL_PYTHON',unset='python3')
hash<-system2(python,c('python_fallback/build_pipeline.py','--manifest-md5'),stdout=TRUE)
stopifnot(length(hash)==1L,grepl('^[a-f0-9]{32}$',hash));release$fallback_bundle_md5<-hash
deploy_json(release,'python_fallback/release.json')
cat('Fallback release bound to unchanged frozen phenotypes. Source hash: ',release$effective_source_hash,'\n',sep='')
