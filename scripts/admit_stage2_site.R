# Coordinator-only: admit a late Stage1 site without changing frozen prototypes.
source('scripts/deployment_lib.R');clif_check_environment()
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=4L)stop('Usage: admit_stage2_site.R EXISTING_BUNDLE STAGE1_RETURN NEW_BUNDLE ACCEPTANCE_JUSTIFICATION')
old<-normalizePath(a[1]);ret<-normalizePath(a[2]);out<-a[3];reason<-a[4]
if(dir.exists(out))stop('Choose a new admission-revision directory.')
bm<-deploy_validate_inventory(old,'bundle_manifest.json');if(bm$source_hash!=deploy_hash()||bm$r_lock_md5!=deploy_lock_hash())stop('Use the frozen release code and environment.')
m<-deploy_read(file.path(ret,'return_manifest.json'));site<-m$site
if(!grepl('^[A-Za-z0-9_-]+$',site)||site%in%unlist(bm$sites)||!nzchar(reason))stop('Require a new site ID and documented acceptance.')
approval<-list(site=site,decision='accepted_legacy_discovery',justification=reason,source_hash=m$source_hash,r_lock_md5=m$r_lock_md5,return_manifest_md5=unname(tools::md5sum(file.path(ret,'return_manifest.json'))))
deploy_stage1_input(ret,approval)
for(f in names(bm$files))deploy_copy(file.path(old,f),file.path(out,f))
for(v in unlist(bm$variants)){
 from<-file.path(old,v,'coordinator/full');to<-file.path(out,v,'coordinator',paste0('loo_',site));deploy_tree(from,to,c('csv','json'))
 fm<-deploy_read(file.path(to,'frozen_manifest.json'));stopifnot(!site%in%unlist(fm$training_sites))
 fm$version<-paste0('loo_',site);fm$held_out_site<-site;fm$late_site_alias_of<-'full';fm$late_site_reason<-'Site did not contribute to original discovery; unchanged full prototypes provide external phenotype mapping.'
 deploy_json(fm,file.path(to,'frozen_manifest.json'))
}
bm$sites<-c(bm$sites,list(site));bm$stage1_acceptance[[site]]<-approval;bm$stage1_manifest_md5[[site]]<-approval$return_manifest_md5
bm$parent_bundle_md5<-unname(tools::md5sum(file.path(old,'bundle_manifest.json')));bm$workspace_id<-basename(out)
bm$late_site_admission<-list(site=site,discovery_contribution=FALSE,all_stage2_analyses=TRUE,full_prototypes_unchanged=TRUE)
bm$files<-deploy_inventory(out);deploy_json(bm,file.path(out,'bundle_manifest.json'));deploy_validate_inventory(out,'bundle_manifest.json')
cat('Admission revision ready: ',out,'\nFull frozen prototypes unchanged. Issue this revision before site runs.\n')
