# Add the new score phase to a frozen site return in a NEW folder.
source('diagnostics/apache_adapted/lib.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=3)stop('Usage: add_to_return.R SITE OLD_SHARE NEW_SHARE')
site<-a[1];old<-a[2];out<-a[3];if(dir.exists(out))stop('Choose fresh destination')
m<-jsonlite::fromJSON(file.path(old,'site_return_manifest.json'),simplifyVector=FALSE);stopifnot(m$site==site)
f<-list.files(old,recursive=TRUE,full.names=FALSE);dir.create(out,recursive=TRUE)
for(k in f){dest<-file.path(out,k);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(file.path(old,k),dest))}
from<-file.path(ap_root(),'sites',site,'export');audit<-jsonlite::fromJSON(file.path(from,'audit.json'));stopifnot(audit$status=='passed')
component<-jsonlite::fromJSON(file.path(from,'component_manifest.json'));stopifnot(component$primary_consensus_id==m$primary_consensus_id)
to<-file.path(out,'federated/apache_adapted');dir.create(to,recursive=TRUE)
for(file in list.files(from,full.names=TRUE)){
 if(grepl('[.]csv$',file)){x<-readr::read_csv(file,show_col_types=FALSE);stopifnot(!any(names(x)%in%c('patient_id','hospitalization_id','recorded_dttm','t0','event_time','sequence')))}
 stopifnot(grepl('[.](csv|json)$',file),file.copy(file,file.path(to,basename(file))))
}
m$parent_return_manifest_md5<-unname(tools::md5sum(file.path(old,'site_return_manifest.json')))
m$apache_adapted_extension_complete<-TRUE;m$apache_not_standard_II_or_IV<-TRUE;m$apache_spec_hash<-fed_hash(ap_spec())
f<-setdiff(list.files(out,recursive=TRUE,full.names=FALSE),'site_return_manifest.json');m$files<-as.list(setNames(unname(tools::md5sum(file.path(out,f))),f))
jsonlite::write_json(m,file.path(out,'site_return_manifest.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
