# One-pass exact aggregate exit histories, including patient-cluster bootstrap.
source('diagnostics/pooled_survival/lib.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=3L)stop('Usage: site_survival.R SITE BUNDLE NEW_OUTPUT')
site<-a[1];bundle<-a[2];out<-a[3];if(dir.exists(out))stop('Choose fresh survival output')
plan<-jsonlite::fromJSON(file.path(bundle,'analysis_plan.json'));B<-plan$prediction$bootstrap
x<-ps_data(site,Sys.getenv('CLIF_FEDERATED_OUTPUT'));patients<-sort(x$patients,method='radix');P<-length(patients)
set.seed(plan$prediction$seed+sum(utf8ToInt(site)));pw<-if(P)replicate(B,tabulate(sample.int(P,P,replace=TRUE),nbins=P)) else matrix(0,0,B)
rows<-list();meta<-list()
for(curve in names(x$sets)){
 d<-x$sets[[curve]]
 status<-ifelse(d$valid_outcome & d$time<=28,d$event_type,'censored')
 stopifnot(all(status%in%c('death','hospice','live_discharge','censored')))
 tick<-ps_tick(pmin(28,d$time));key<-paste(tick,status,sep='|');keys<-unique(key)
 h<-tibble::tibble(time_tick=as.numeric(sub('[|].*','',keys)),exit_type=sub('.*[|]','',keys),n=tabulate(match(key,keys),nbins=length(keys)))
 for(b in seq_len(B)){
  w<-pw[match(d$patient_id,patients),b];counts<-numeric(length(keys))
  if(length(keys)){q<-rowsum(w,match(key,keys),reorder=FALSE);counts[as.integer(rownames(q))]<-q[,1]}
  h[[sprintf('n_b%03d',b)]]<-counts
 }
 h$site<-rep(site,nrow(h));h$curve<-rep(curve,nrow(h));rows[[curve]]<-h
 meta[[curve]]<-tibble::tibble(site=site,curve=curve,n=nrow(d),patients=dplyr::n_distinct(d$patient_id),unknown_outcomes=sum(!d$valid_outcome),max_followup=if(nrow(d))min(28,max(d$time)) else 0)
}
ps_export(dplyr::bind_rows(rows),file.path(out,'exit_histograms.csv'));ps_export(dplyr::bind_rows(meta),file.path(out,'cohort_counts.csv'))
f<-list.files(out,full.names=TRUE)
jsonlite::write_json(list(status='passed',site=site,consensus_id=jsonlite::fromJSON(file.path(Sys.getenv('CLIF_FEDERATED_OUTPUT'),'coordinator/full/frozen_manifest.json'))$consensus_id,bootstrap=B,method='Exact relative-duration exit counts; event and censor histories permit central reconstruction at the union of site times without a second site pass',bootstrap_unit='patients within site; all admissions retained',unknown_outcomes='censored at observed exit',files=as.list(setNames(unname(tools::md5sum(f)),basename(f)))),file.path(out,'manifest.json'),auto_unbox=TRUE,pretty=TRUE)
