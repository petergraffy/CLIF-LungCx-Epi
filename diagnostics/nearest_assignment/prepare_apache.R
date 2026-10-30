source('diagnostics/apache_adapted/lib.R')
parent<-Sys.getenv('CLIF_NEAREST_APACHE_PARENT',unset='output/run_APACHE_adapted_20261030_v2')
if(dir.exists(ap_root()))stop('Choose fresh adapted-score output')
consensus<-jsonlite::fromJSON(file.path(ap_primary(),'coordinator/full/frozen_manifest.json'))$consensus_id
for(site in ap_sites()){
 old<-file.path(parent,'sites',site);dest<-file.path(ap_root(),'sites',site)
 dir.create(file.path(dest,'local'),recursive=TRUE);dir.create(file.path(dest,'export'))
 for(f in c('components.rds','source_metadata_match_audit.csv'))stopifnot(file.copy(file.path(old,'local',f),file.path(dest,'local',f)))
 a<-readRDS(file.path(dest,'local/components.rds'));d<-readRDS(file.path(ap_primary(),'sites',site,'local/full/outcomes.rds'))$data
 stopifnot(identical(a$hospitalization_id,d$hospitalization_id))
 rows<-lapply(c(APACHE_COMPONENTS,'adapted_total','observed_hct_total'),function(k){v<-a[[paste0('apache_',k)]];n<-sum(is.finite(v));miss<-sum(!is.finite(v));tibble::tibble(site=site,component=k,n=nrow(a),observed_n=n,missing_n=miss,observed_lower=n,observed_upper=n,missing_lower=miss,missing_upper=miss)})
 ap_export(dplyr::bind_rows(rows),file.path(dest,'export/component_completeness.csv'))
 flags<-c('hb_proxy_used','map_derived','gcs_derived','bicarbonate_substituted')
 ap_export(dplyr::bind_rows(lapply(flags,function(k){n<-sum(a[[paste0('apache_',k)]]);tibble::tibble(site=site,flag=k,n=n,n_lower=n,n_upper=n)})),file.path(dest,'export/provenance_counts.csv'))
 patterns<-tibble::tibble(pattern=apply(is.na(a[paste0('apache_',APACHE_COMPONENTS)]),1,paste0,collapse='')) |> dplyr::count(pattern,name='n')
 patterns$n_lower<-patterns$n_upper<-patterns$n;patterns$site<-site;ap_export(patterns,file.path(dest,'export/missing_patterns.csv'))
 m<-jsonlite::fromJSON(file.path(old,'export/component_manifest.json'),simplifyVector=FALSE)
 m$parent_primary_consensus_id<-m$primary_consensus_id;m$primary_consensus_id<-consensus;m$component_values_reused_unchanged<-TRUE;m$specification<-ap_spec();m$spec_hash<-fed_hash(ap_spec());m$parent_component_manifest_md5<-unname(tools::md5sum(file.path(old,'export/component_manifest.json')))
 fed_json(m,file.path(dest,'export/component_manifest.json'))
}
cat('Adapted-score component values reused; phenotype-dependent models will be refit\n')
