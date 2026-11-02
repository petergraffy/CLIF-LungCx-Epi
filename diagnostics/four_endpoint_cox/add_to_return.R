# Append the aggregate four-endpoint package to an immutable prior site return.
source('code/lib/load.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=4)stop('Usage: add_to_return.R SITE FOUR_COX_SHARE OLD_RETURN NEW_RETURN')
site<-a[1];extension<-a[2];old<-a[3];out<-a[4]
if(dir.exists(out))stop('Choose fresh output')
m<-jsonlite::fromJSON(file.path(old,'site_return_manifest.json'),simplifyVector=FALSE)
e<-jsonlite::fromJSON(file.path(extension,'manifest.json'))
stopifnot(m$site==site,e$site==site,e$status=='passed',e$primary_consensus_id==m$primary_consensus_id,setequal(e$endpoints,cox_endpoints()),e$original_endpoints_reproduced,e$assignments_unchanged)
for(k in names(e$files)){
 stopifnot(!grepl('(^/|(^|/)\\.\\.(/|$))',k),unname(tools::md5sum(file.path(extension,k)))==e$files[[k]])
}
dir.create(out,recursive=TRUE)
for(k in list.files(old,recursive=TRUE,full.names=FALSE)){
 p<-file.path(out,k);dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(file.path(old,k),p))
}
for(k in list.files(extension,recursive=TRUE,full.names=FALSE)){
 p<-file.path(out,'federated/four_endpoint_cox',k);dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(file.path(extension,k),p))
}
f<-setdiff(list.files(out,recursive=TRUE,full.names=FALSE),'site_return_manifest.json')
readr::write_csv(tibble::tibble(site=site,file=f,status='captured'),file.path(out,'supplement/return_file_inventory.csv'))
m$parent_return_manifest_md5<-unname(tools::md5sum(file.path(old,'site_return_manifest.json')))
m$four_endpoint_cox_complete<-TRUE;m$four_endpoint_cox_endpoints<-cox_endpoints()
m$four_endpoint_cox_primary_export<-'federated/four_endpoint_cox/primary/'
m$files<-as.list(setNames(unname(tools::md5sum(file.path(out,f))),f))
jsonlite::write_json(m,file.path(out,'site_return_manifest.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
cat(site,'four-endpoint Cox primary/sensitivity returns bundled\n')
