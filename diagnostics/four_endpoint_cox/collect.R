# Coordinator: aggregate-only association meta-analysis, never clinical caches.
source('code/lib/load.R');source('code/lib/central.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)<2)stop('Usage: collect.R NEW_OUTPUT SITE_SHARE [...]')
out<-a[1];shares<-a[-1];if(dir.exists(out))stop('Choose fresh central output')
m<-lapply(shares,function(p)jsonlite::fromJSON(file.path(p,'manifest.json')))
sites<-vapply(m,function(x)x$site,character(1));stopifnot(!anyDuplicated(sites),all(vapply(m,function(x)x$status=='passed'&&x$original_endpoints_reproduced&&x$assignments_unchanged,logical(1))),length(unique(vapply(m,function(x)x$primary_consensus_id,character(1))))==1)
variants<-m[[1]]$variants;stopifnot(all(vapply(m,function(x)identical(x$variants,variants)&&setequal(x$endpoints,cox_endpoints()),logical(1))))
for(i in seq_along(shares))for(f in names(m[[i]]$files))stopifnot(unname(tools::md5sum(file.path(shares[i],f)))==m[[i]]$files[[f]])
dir.create(out,recursive=TRUE);all_assoc<-all_meta<-all_audit<-list()
for(variant in variants){
 stopifnot(length(unique(vapply(m,function(x)x$inputs[[variant]]$consensus_id,character(1))))==1)
 for(name in c('associations','model_metadata','association_covariances','event_partitions','model_PH_and_denominators','model_PH_terms','exploratory_piecewise_associations','exploratory_piecewise_metadata')){
  rows<-lapply(seq_along(shares),function(i){
   z<-readr::read_csv(file.path(shares[i],variant,paste0(name,'.csv')),show_col_types=FALSE)
   stopifnot(!any(names(z)%in%c('patient_id','hospitalization_id','t0','event_time','sequence')))
   if(nrow(z))stopifnot(all(z$site==sites[i]),all(z$analysis_variant==variant))
   z
  })
  z<-dplyr::bind_rows(rows);dir.create(file.path(out,variant),recursive=TRUE,showWarnings=FALSE)
  readr::write_csv(z,file.path(out,variant,paste0(name,'_by_site.csv')),na='')
  if(name=='associations')assoc<-z
  if(name=='exploratory_piecewise_associations')piece_assoc<-z
  if(name=='model_metadata'){
   stopifnot(all(vapply(split(z,z$site),function(q)setequal(q$endpoint,cox_endpoints())&&nrow(q)==12,logical(1))))
   all_meta[[variant]]<-z
  }
 }
 pool<-assoc |> dplyr::group_by(model,endpoint,phenotype,reference,effect_type) |> dplyr::group_modify(~pool_effect(.x)) |> dplyr::ungroup()
 pool$analysis_variant<-variant
 readr::write_csv(pool,file.path(out,variant,'pooled_associations.csv'),na='')
 piece_pool<-piece_assoc |> dplyr::group_by(interval_start_day,interval_end_day,endpoint,phenotype,reference,analysis) |> dplyr::group_modify(~pool_effect(.x)) |> dplyr::ungroup()
 piece_pool$analysis_variant<-variant
 readr::write_csv(piece_pool,file.path(out,variant,'pooled_exploratory_piecewise_associations.csv'),na='')
 all_assoc[[variant]]<-pool
}
readr::write_csv(dplyr::bind_rows(all_assoc),file.path(out,'pooled_associations_all_variants.csv'),na='')
readr::write_csv(dplyr::bind_rows(all_meta),file.path(out,'model_metadata_all_variants.csv'),na='')
primary<-all_assoc$primary |> dplyr::filter(model=='core')
readr::write_csv(primary,file.path(out,'table2_four_endpoint_core_cox.csv'),na='')
jsonlite::write_json(list(status='passed',sites=sites,primary_consensus_id=m[[1]]$primary_consensus_id,variants=variants,endpoints=cox_endpoints(),patient_level_inputs_read=FALSE,method='Site cause-specific Cox with patient-clustered SE; REML modified Hartung-Knapp pooling; >=2 estimable sites required',horizon='28 days after each analysis landmark',interpretation='Live-discharge HR >1 is faster live non-hospice discharge among those still hospitalized; not an absolute probability or inverse mortality HR',source_manifest_md5=as.list(setNames(unname(tools::md5sum(file.path(shares,'manifest.json'))),sites))),file.path(out,'manifest.json'),auto_unbox=TRUE,pretty=TRUE)
cat('Four-endpoint primary and sensitivity Cox meta-analysis passed\n')
