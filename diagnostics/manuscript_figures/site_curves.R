# One-site worker: only this site's frozen clinical cache is opened.
source('code/lib/load.R')
source('diagnostics/manuscript_figures/cif_lib.R')
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=3)stop('Usage: site_curves.R SITE PRIMARY_ROOT NEW_SITE_OUTPUT')
site<-a[1];primary<-a[2];out<-a[3];if(dir.exists(out))stop('Choose fresh site output')
cache<-file.path(primary,'sites',site,'local/full/outcomes.rds')
d<-readRDS(cache)$data
consensus<-jsonlite::fromJSON(file.path(primary,'coordinator/full/frozen_manifest.json'))$consensus_id
cif<-released_cif_grid(d,0:28,10L) |> dplyr::mutate(site=site,version='full',consensus_id=consensus)
stopifnot(!any(names(cif)%in%c('patient_id','hospitalization_id','event_time','t0','time_from_landmark')))
for(ph in unique(cif$phenotype))for(day in unique(cif$day[cif$phenotype==ph])){
 z<-cif[cif$phenotype==ph&cif$day==day,]
 if(all(z$cell_status=='observed'))stopifnot(abs(sum(z$estimate)-1)<1e-8)
}
# Check all previously released checkpoints are exactly reproduced.
old<-readr::read_csv(file.path(primary,'sites',site,'export/outcomes/full/cumulative_incidence.csv'),show_col_types=FALSE)
z<-dplyr::inner_join(old,cif,by=c('phenotype','day','event'),suffix=c('_old','_new'))
stopifnot(nrow(z)==nrow(old),all(is.finite(z$estimate_new[is.finite(z$estimate_old)])),all(abs(z$estimate_old-z$estimate_new)<1e-12,na.rm=TRUE))
dir.create(out,recursive=TRUE);dest<-file.path(out,'share');dir.create(dest)
readr::write_csv(cif,file.path(dest,'cumulative_incidence_daily.csv'),na='')
jsonlite::write_json(list(status='passed',site=site,consensus_id=consensus,grid_days=0:28,
 method='Unadjusted Aalen-Johansen; death, hospice, live discharge mutually exclusive; hospitalization clock from72h; patient-clustered SE',
 aggregate_policy='exact counts and estimates; no disclosure suppression',
 existing_checkpoint_estimates_unchanged=TRUE,no_patient_identifiers_exported=TRUE,
 input_cache_md5=unname(tools::md5sum(cache)),files=as.list(setNames(unname(tools::md5sum(file.path(dest,'cumulative_incidence_daily.csv'))),'cumulative_incidence_daily.csv'))),
 file.path(dest,'curve_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
cat(site,'daily Aalen-Johansen export passed; frozen checkpoints reproduced\n')
