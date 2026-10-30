# Version a diagnostic extension without overwriting an existing site return.
source('code/lib/load.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=3)stop('Usage: add_model_diagnostics.R RETURN_SHARE DIAGNOSTIC_SHARE NEW_SHARE')
old<-a[1];extension<-a[2];out<-a[3]
if(dir.exists(out))stop('Choose a fresh destination')
m<-jsonlite::fromJSON(file.path(old,'site_return_manifest.json'),simplifyVector=FALSE)
e<-jsonlite::fromJSON(file.path(extension,'manifest.json'))
stopifnot(e$status=='passed',e$site==m$site,isTRUE(e$original_coefficients_reproduced))
f<-list.files(old,recursive=TRUE,full.names=FALSE);dir.create(out,recursive=TRUE)
for(k in f){dest<-file.path(out,k);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(file.path(old,k),dest))}
for(k in c('model_PH_and_denominators.csv','model_PH_terms.csv','exploratory_piecewise_associations.csv','exploratory_piecewise_metadata.csv')){
 x<-readr::read_csv(file.path(extension,k),show_col_types=FALSE)
 stopifnot(all(x$site==e$site),!any(names(x)%in%c('patient_id','hospitalization_id','t0','event_time','sequence')))
 readr::write_csv(x,file.path(out,'federated',k),na='')
}
m$parent_return_manifest_md5<-unname(tools::md5sum(file.path(old,'site_return_manifest.json')))
m$diagnostic_extension_manifest_md5<-unname(tools::md5sum(file.path(extension,'manifest.json')))
m$diagnostic_extension<-list(original_coefficients_reproduced=TRUE,original_models_not_modified=TRUE,diagnostic_refits_local_only=TRUE)
f<-setdiff(list.files(out,recursive=TRUE,full.names=FALSE),'site_return_manifest.json')
m$files<-as.list(setNames(unname(tools::md5sum(file.path(out,f))),f))
jsonlite::write_json(m,file.path(out,'site_return_manifest.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
