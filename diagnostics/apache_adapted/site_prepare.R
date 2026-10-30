source('diagnostics/apache_adapted/lib.R');source('diagnostics/sofa_completeness/context.R')
site<-commandArgs(trailingOnly=TRUE)[1];stopifnot(site%in%ap_sites());home<-file.path(ap_root(),'sites',site)
dir.create(file.path(home,'local'),recursive=TRUE,showWarnings=FALSE);dir.create(file.path(home,'export'),recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(home,'local/components.rds')))stop('Components already exist')
dis<-readRDS(file.path(ap_primary(),'sites',site,'local/discovery.rds'));m<-jsonlite::fromJSON(file.path(dis$input_root,'local/input_manifest.json'))
p<-readRDS(file.path(m$baseline_root,'local/prepared.rds'));d<-readRDS(file.path(ap_primary(),'sites',site,'local/full/outcomes.rds'))$data
source_path<-Sys.getenv(paste0('CLIF_APACHE_TABLES_',site),unset=dirname(m$source_file))
raw<-sofa_source_rows(source_path,'labs',d$hospitalization_id)
cached<-p$tables$labs[p$tables$labs$hospitalization_id%in%d$hospitalization_id,]
enriched<-apache_enrich_cached_labs(cached,raw);p$tables$labs<-enriched$labs
readr::write_csv(enriched$audit,file.path(home,'local/source_metadata_match_audit.csv'))
a<-apache_adapted_features(d,p$tables);stopifnot(identical(a$hospitalization_id,d$hospitalization_id),all(a$apache_adapted_complete==!is.na(a$apache_adapted_total)))
saveRDS(a,file.path(home,'local/components.rds'))
rows<-list()
for(k in c(APACHE_COMPONENTS,'adapted_total','observed_hct_total')){
 v<-a[[paste0('apache_',k)]];b<-fed_band_n(sum(is.finite(v)));miss<-fed_band_n(sum(!is.finite(v)))
 rows[[k]]<-tibble::tibble(site=site,component=k,n=nrow(a),observed_lower=b[1],observed_upper=b[2],missing_lower=miss[1],missing_upper=miss[2])
}
ap_export(dplyr::bind_rows(rows),file.path(home,'export/component_completeness.csv'))
flags<-c('hb_proxy_used','map_derived','gcs_derived','bicarbonate_substituted')
flags<-dplyr::bind_rows(lapply(flags,function(k){b<-fed_band_n(sum(a[[paste0('apache_',k)]]));tibble::tibble(site=site,flag=k,n_lower=b[1],n_upper=b[2])}))
ap_export(flags,file.path(home,'export/provenance_counts.csv'))
patterns<-a |> dplyr::mutate(pattern=apply(is.na(a[paste0('apache_',APACHE_COMPONENTS)]),1,paste0,collapse='')) |> dplyr::count(pattern,name='n')
b<-t(vapply(patterns$n,fed_band_n,numeric(2)));patterns$n_lower<-b[,1];patterns$n_upper<-b[,2];patterns$n<-NULL;patterns$site<-site
ap_export(patterns,file.path(home,'export/missing_patterns.csv'))
fed_json(list(site=site,n=nrow(a),specification=ap_spec(),spec_hash=fed_hash(ap_spec()),source_id=p$source_id,score_source_md5=as.list(tools::md5sum('code/lib/apache_adapted.R')),primary_consensus_id=jsonlite::fromJSON(file.path(ap_primary(),'coordinator/full/frozen_manifest.json'))$consensus_id,not_standard_APACHE_II=TRUE,not_APACHE_IV=TRUE),file.path(home,'export/component_manifest.json'))
cat(site,'components complete; observed adapted total',sum(a$apache_adapted_complete),'of',nrow(a),'\n')
