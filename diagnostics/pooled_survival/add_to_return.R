# Bundle the two local aggregate passes; never copy a local clinical cache.
source('code/lib/load.R')
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=5)stop('Usage: add_to_return.R SITE TIME_SHARE RISK_SHARE OLD_RETURN NEW_RETURN')
site<-a[1];times<-a[2];risks<-a[3];old<-a[4];out<-a[5]
if(dir.exists(out))stop('Choose fresh return output')
m<-jsonlite::fromJSON(file.path(old,'site_return_manifest.json'),simplifyVector=FALSE)
tm<-jsonlite::fromJSON(file.path(times,'manifest.json'));rm<-jsonlite::fromJSON(file.path(risks,'manifest.json'))
stopifnot(m$site==site,tm$site==site,rm$site==site,tm$status=='passed',rm$status=='passed',tm$consensus_id==m$primary_consensus_id,rm$consensus_id==m$primary_consensus_id,rm$bootstrap==200,rm$local_AJ_reconstruction_passed,rm$local_KM_reconstruction_passed)
for(p in c(times,risks)){
 em<-jsonlite::fromJSON(file.path(p,'manifest.json'))
 for(k in names(em$files))stopifnot(unname(tools::md5sum(file.path(p,k)))==em$files[[k]])
}
cohorts<-readr::read_csv(file.path(times,'cohort_counts.csv'),show_col_types=FALSE)
stopifnot(all(c('death','death_or_hospice','live_discharge','any_exit')%in%cohorts$endpoint))
dir.create(out,recursive=TRUE)
for(f in list.files(old,recursive=TRUE,full.names=FALSE)){
 dest<-file.path(out,f);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
 stopifnot(file.copy(file.path(old,f),dest))
}
for(phase in c('times','risksets')){
 p<-if(phase=='times')times else risks
 dest<-file.path(out,'federated/pooled_survival',phase);dir.create(dest,recursive=TRUE)
 for(f in list.files(p,full.names=TRUE))stopifnot(file.copy(f,file.path(dest,basename(f))))
}
f<-setdiff(list.files(out,recursive=TRUE,full.names=FALSE),'site_return_manifest.json')
readr::write_csv(tibble::tibble(site=site,file=f,status='captured'),file.path(out,'supplement/return_file_inventory.csv'))
m$parent_return_manifest_md5<-unname(tools::md5sum(file.path(old,'site_return_manifest.json')))
m$pooled_survival_discharge_complete<-TRUE;m$pooled_grid_md5<-rm$grid_md5
m$files<-as.list(setNames(unname(tools::md5sum(file.path(out,f))),f))
jsonlite::write_json(m,file.path(out,'site_return_manifest.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
cat(site,'pooled survival/live-discharge aggregate returns bundled and hashes verified\n')
