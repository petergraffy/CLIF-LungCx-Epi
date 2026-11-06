# Validate native and fallback packets against their actual, distinct contracts.
verify_stage2_return<-function(m,bm,bhash){
 release<-deploy_read('python_fallback/release.json')
 stopifnot(m$stage=='stage2',m$status=='complete',m$workflow==bm$workflow,m$primary_consensus_id==bm$primary_consensus_id,m$r_lock_md5==bm$r_lock_md5)
 if(is.null(m$python_fallback)){
  stopifnot(m$parameter_bundle_md5==bhash,m$source_hash==bm$source_hash)
 }else{
  p<-m$python_fallback
  stopifnot(bhash==release$parent_bundle_md5,bm$source_hash==release$parent_source_hash,
   p$release==release$release,p$parquet_backend=='pyarrow',identical(p$r_arrow_required,FALSE),
   p$parent_bundle_md5==bhash,p$parent_source_hash==bm$source_hash,
   m$parameter_bundle_md5==release$fallback_bundle_md5,m$source_hash==release$effective_source_hash,
   setequal(names(m$source_files_md5),names(release$effective_source_files_md5)),
   identical(m$source_files_md5[names(release$effective_source_files_md5)],release$effective_source_files_md5))
 }
 invisible(TRUE)
}
