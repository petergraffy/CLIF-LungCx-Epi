source('code/lib/load.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=4)stop('Usage: add_survival_to_return.R SITE SURVIVAL_SHARE OLD_SHARE NEW_SHARE')
site<-a[1];extension<-a[2];old<-a[3];out<-a[4];if(dir.exists(out))stop('Choose fresh return')
m<-jsonlite::fromJSON(file.path(old,'site_return_manifest.json'),simplifyVector=FALSE)
em<-jsonlite::fromJSON(file.path(extension,'manifest.json'))
stopifnot(m$site==site,em$site==site,em$status=='passed',em$primary_consensus_id==m$primary_consensus_id)
for(k in names(em$files))stopifnot(unname(tools::md5sum(file.path(extension,k)))==em$files[[k]])
dir.create(out,recursive=TRUE)
for(f in list.files(old,recursive=TRUE,full.names=FALSE)){dest<-file.path(out,f);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(file.path(old,f),dest))}
dest<-file.path(out,'federated/survival_tdauc');dir.create(dest)
for(f in list.files(extension,full.names=TRUE))stopifnot(file.copy(f,file.path(dest,basename(f))))
m$parent_return_manifest_md5<-unname(tools::md5sum(file.path(old,'site_return_manifest.json')));m$survival_tdauc_complete<-TRUE
f<-setdiff(list.files(out,recursive=TRUE,full.names=FALSE),'site_return_manifest.json');m$files<-as.list(setNames(unname(tools::md5sum(file.path(out,f))),f))
jsonlite::write_json(m,file.path(out,'site_return_manifest.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
