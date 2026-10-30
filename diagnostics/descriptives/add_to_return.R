source('code/lib/load.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=5)stop('Usage: add_to_return.R SITE DESCRIPTORS_ROOT BASELINE_ROOT OLD_SHARE NEW_SHARE')
site<-a[1];desc<-file.path(a[2],'share/descriptives');baseline<-a[3];old<-a[4];out<-a[5]
if(dir.exists(out))stop('Choose fresh return folder')
m<-jsonlite::fromJSON(file.path(old,'site_return_manifest.json'),simplifyVector=FALSE);d<-jsonlite::fromJSON(file.path(desc,'manifest.json'))
stopifnot(m$site==site,d$site==site,d$status=='passed',d$grouping=='phenotype',d$consensus_id==m$primary_consensus_id)
dir.create(out,recursive=TRUE)
for(f in list.files(old,recursive=TRUE,full.names=FALSE)){dest<-file.path(out,f);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(file.path(old,f),dest))}
to<-file.path(out,'federated/descriptives');dir.create(to,recursive=TRUE)
for(f in list.files(desc,full.names=TRUE))stopifnot(file.copy(f,file.path(to,basename(f))))
to<-file.path(out,'federated/cohort_audit');dir.create(to,recursive=TRUE)
for(f in c('cohort_exclusions.csv','cohort_by_year.csv','cohort_patient_counts.csv','outcome_timestamp_audit.csv','source_availability.csv','analysis_plan.csv','arf_detection_audit.csv','derivation_status.csv')){
 src<-file.path(baseline,'share',f);stopifnot(file.exists(src),file.copy(src,file.path(to,f)))
}
m$parent_return_manifest_md5<-unname(tools::md5sum(file.path(old,'site_return_manifest.json')));m$publication_descriptives_complete<-TRUE;m$publication_table2_grouping<-'frozen_shared_primary'
f<-setdiff(list.files(out,recursive=TRUE,full.names=FALSE),'site_return_manifest.json');m$files<-as.list(setNames(unname(tools::md5sum(file.path(out,f))),f))
jsonlite::write_json(m,file.path(out,'site_return_manifest.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
