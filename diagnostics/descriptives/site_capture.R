source('diagnostics/apache_adapted/lib.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=2)stop('Usage: site_capture.R SITE NEW_OUTPUT')
site<-a[1];out<-a[2];if(dir.exists(out))stop('Choose a fresh descriptor root')
home<-file.path(ap_primary(),'sites',site);disc<-readRDS(file.path(home,'local/discovery.rds'))
m<-jsonlite::fromJSON(file.path(disc$input_root,'local/input_manifest.json'));p<-readRDS(file.path(m$baseline_root,'local/prepared.rds'));d<-readRDS(file.path(home,'local/full/outcomes.rds'))$data
cfg<-list(site_name=site,analysis=p$parameters)
apfile<-file.path(ap_root(),'sites',site,'local/components.rds');ap<-if(file.exists(apfile))readRDS(apfile) else NULL
phenotypes<-rownames(readRDS(file.path(ap_primary(),'coordinator/full/prototypes.rds')))
consensus<-jsonlite::fromJSON(file.path(ap_primary(),'coordinator/full/frozen_manifest.json'))$consensus_id
capture_descriptive_tables(p,d,cfg,out,disc$panel,phenotypes,ap,consensus)
cat(site,'Table1/Table2 descriptor capture passed\n')
