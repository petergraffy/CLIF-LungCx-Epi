# Reuse the learned templates/local clusters; version only the assignment rule.
source('diagnostics/federated/contract.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)<2)stop('Usage: prepare.R PARENT_ROOT NEW_ROOT [VARIANT]')
parent<-normalizePath(a[1]);out<-a[2];if(dir.exists(out))stop('Choose fresh nearest-assignment root')
if(length(a)==3){
 Sys.setenv(CLIF_REVISION_VARIANT=a[3])
 source(if(a[3]=='absorbing_v2')'diagnostics/revision20261022/absorbing_contract.R' else 'diagnostics/revision20261022/variant_contract.R')
}
dir.create(out,recursive=TRUE);parents<-list()
versions<-intersect(c('full','loo_UCMC','loo_NU','loo_RUSH'),list.dirs(file.path(parent,'coordinator'),recursive=FALSE,full.names=FALSE))
stopifnot('full'%in%versions)
for(version in versions){
 old<-file.path(parent,'coordinator',version);dest<-file.path(out,'coordinator',version);dir.create(dest,recursive=TRUE)
 for(f in c('prototypes.rds','prototypes.csv','candidate_solutions.csv','local_cluster_diagnostic_mapping.csv'))if(file.exists(file.path(old,f)))stopifnot(file.copy(file.path(old,f),file.path(dest,f)))
 m<-jsonlite::fromJSON(file.path(old,'frozen_manifest.json'),simplifyVector=FALSE);proto<-readRDS(file.path(dest,'prototypes.rds'))
 m$parent_consensus_id<-m$consensus_id;m$parent_frozen_manifest_md5<-unname(tools::md5sum(file.path(old,'frozen_manifest.json')))
 m$template_learning_specification<-m$specification;m$schema<-FED_SCHEMA;m$specification<-fed_spec();m$spec_hash<-fed_hash(fed_spec())
 m$consensus_id<-substr(fed_hash(list(spec=fed_spec(),training_sites=as.character(unlist(m$training_sites)),source=m$source_hashes,prototypes=proto)),1,20)
 m$assignment<-'Always closest weighted profile distance for information-sufficient stays; exact ties use frozen template order; distance/margin diagnostic only'
 m$template_matrix_unchanged<-TRUE;m$template_parent<-old
 m$assignment_source_md5<-as.list(tools::md5sum('diagnostics/federated/contract.R'))
 jsonlite::write_json(m,file.path(dest,'frozen_manifest.json'),auto_unbox=TRUE,pretty=TRUE,null='null',digits=16)
 parents[[version]]<-list(parent=m$parent_consensus_id,current=m$consensus_id,prototype_md5=unname(tools::md5sum(file.path(dest,'prototypes.rds'))))
}
for(site in c('UCMC','NU','RUSH')){
 dest<-file.path(out,'sites',site,'local');dir.create(dest,recursive=TRUE)
 for(f in c('discovery.rds','patient_features.rds','discovery_audit.json')){
  from<-file.path(parent,'sites',site,'local',f);if(file.exists(from))stopifnot(file.copy(from,file.path(dest,f)))
 }
 # Current aggregate views are exact, while original source hashes above retain
 # the immutable coarsened artifacts from which the templates were learned.
 old<-file.path(parent,'sites',site,'export/discovery');ex<-file.path(out,'sites',site,'export/discovery');dir.create(ex,recursive=TRUE)
 for(f in list.files(old,full.names=TRUE))stopifnot(file.copy(f,file.path(ex,basename(f))))
 d<-readRDS(file.path(dest,'discovery.rds'))
 counts<-d$data |> dplyr::count(local_cluster,name='exact_n')
 z<-readr::read_csv(file.path(ex,'clusters.csv'),show_col_types=FALSE) |> dplyr::left_join(counts,by='local_cluster')
 z$n<-z$n_lower<-z$n_upper<-z$exact_n;z$exact_n<-NULL
 fed_write(z,file.path(ex,'clusters.csv'))
 f<-readr::read_csv(file.path(ex,'fingerprints.csv'),show_col_types=FALSE)
 modal<-fed_modal(d$panel) |> dplyr::left_join(d$data |> dplyr::select(hospitalization_id,local_cluster),by='hospitalization_id') |> dplyr::count(local_cluster,domain,bin,state,name='count')
 f<-f |> dplyr::select(-dplyr::any_of(c('count','denominator'))) |> dplyr::left_join(modal,by=c('local_cluster','domain','bin','state')) |> dplyr::left_join(counts,by='local_cluster')
 f$count[is.na(f$count)]<-0;f$denominator<-f$exact_n;f$exact_n<-NULL;fed_write(f,file.path(ex,'fingerprints.csv'))
 dm<-jsonlite::fromJSON(file.path(ex,'manifest.json'),simplifyVector=FALSE);dm$parent_discovery_root<-old;dm$reused_local_clustering<-TRUE;dm$aggregate_policy<-'exact aggregate counts; original template intervals retained';dm$assignment_schema<-FED_SCHEMA
 jsonlite::write_json(dm,file.path(ex,'manifest.json'),auto_unbox=TRUE,pretty=TRUE,null='null',digits=16)
}
fed_json(list(status='prepared',parent_root=parent,template_learning_rerun=FALSE,local_clustering_rerun=FALSE,assignment_rule=fed_spec()$assignment_rule,versions=parents),file.path(out,'nearest_assignment_parent.json'))
cat('Frozen local clusters/templates retained; new nearest-assignment contracts prepared\n')
