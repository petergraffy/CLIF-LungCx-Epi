# Second local pass: evaluate this site's risk sets on the common public grid.
source('diagnostics/pooled_survival/lib.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=4)stop('Usage: site_risksets.R SITE PRIMARY_ROOT CENTRAL_GRID NEW_OUTPUT')
site<-a[1];primary<-a[2];center<-a[3];out<-a[4];if(dir.exists(out))stop('Choose fresh site risk-set output')
cm<-jsonlite::fromJSON(file.path(center,'manifest.json'));stopifnot(site%in%cm$sites,unname(tools::md5sum(file.path(center,'shared_grid.csv')))==cm$grid_md5)
grid<-readr::read_csv(file.path(center,'shared_grid.csv'),show_col_types=FALSE)
x<-ps_data(site,primary);keys<-unique(grid[c('curve','endpoint')]);nominal<-list();cache<-list()
for(i in seq_len(nrow(keys))){
 key<-keys[i,];g<-grid[grid$curve==key$curve&grid$endpoint==key$endpoint,];d<-x$sets[[key$curve]]
 event<-ps_event(d,key$endpoint);q<-ps_counts(d$time,event,g$time)
 # This local reconstruction must equal a conventional local KM fit.
 if(nrow(d)){
  fit<-survival::survfit(survival::Surv(d$time,event)~1,timefix=FALSE)
  expected<-summary(fit,times=ps_time(g$time_tick),extend=TRUE)$surv
  stopifnot(all(abs(ps_product(q$risk,q$events)-expected)<1e-10))
 }
 nominal[[i]]<-g |> dplyr::mutate(site=site,n_risk=q$risk,n_event=q$events,n=nrow(d),max_followup=if(nrow(d))max(d$time) else 0)
 cache[[i]]<-list(g=g,d=d,event=event)
}
# Independent local AJ reconstruction before exporting any aggregates.
for(curve_key in names(x$sets)){
 d<-x$sets[[curve_key]];if(!nrow(d))next
 g<-grid[grid$curve==curve_key&grid$endpoint=='live_discharge',]
 live<-ps_counts(d$time,ps_event(d,'live_discharge'),g$time)
 exits<-ps_counts(d$time,ps_event(d,'any_exit'),g$time)
 status<-ifelse(d$valid_outcome,d$event_type,'censored')
 fit<-survival::survfit(survival::Surv(d$time,factor(status,levels=c('censored','death','hospice','live_discharge')))~1,timefix=FALSE)
 expected<-summary(fit,times=ps_time(g$time_tick),extend=TRUE)$pstate[,'live_discharge']
 stopifnot(all(abs(ps_incidence(live$risk,live$events,exits$events)-expected)<1e-10))
}
nom<-dplyr::bind_rows(nominal);ps_export(nom,file.path(out,'risksets.csv'))
B<-200L;patients<-x$patients;P<-length(patients);weights<-matrix(0,P,B)
set.seed(20261030+sum(utf8ToInt(site)))
for(b in seq_len(B))weights[,b]<-tabulate(sample.int(P,P,replace=TRUE),nbins=P)
boots<-list()
for(i in seq_along(cache)){
 z<-cache[[i]];d<-z$d;idx<-match(d$patient_id,patients)
 for(b in seq_len(B)){
  w<-weights[idx,b];q<-ps_counts(d$time,z$event,z$g$time,w)
  maxfu<-if(any(w>0))max(d$time[w>0]) else 0
  boots[[paste(i,b)]]<-tibble::tibble(site=site,grid_id=z$g$grid_id,replicate=b,n_risk=q$risk,n_event=q$events,max_followup=maxfu)
 }
}
ps_export(dplyr::bind_rows(boots),file.path(out,'bootstrap_risksets.csv'))
files<-list.files(out,full.names=TRUE)
jsonlite::write_json(list(status='passed',site=site,consensus_id=cm$consensus_id,grid_md5=cm$grid_md5,bootstrap=200,bootstrap_unit='patients within site; all their admissions retained; same weights across curves/endpoints',local_KM_reconstruction_passed=TRUE,local_AJ_reconstruction_passed=TRUE,clinical_models_refitted=FALSE,
 files=as.list(setNames(unname(tools::md5sum(files)),basename(files)))),file.path(out,'manifest.json'),auto_unbox=TRUE,pretty=TRUE)
cat(site,'risk sets and200 patient-cluster bootstrap aggregates exported\n')
