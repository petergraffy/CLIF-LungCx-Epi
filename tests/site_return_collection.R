source('code/lib/load.R')
base<-tempfile('site_collection_');dir.create(base)
required<-c('supplement/return_file_inventory.csv','supplement/sofa_original_components.csv','supplement/sofa_recovery_components.csv','supplement/sofa_source_audit_status.csv','supplement/site_reporting_metadata.csv','federated/primary/discovery/silhouette.csv','federated/primary/outcomes/full/model_metadata.csv','federated/scores/transport_metrics.csv','federated/absorbing/early_events.csv','federated/model_PH_and_denominators.csv','federated/model_PH_terms.csv','federated/exploratory_piecewise_associations.csv','federated/exploratory_piecewise_metadata.csv',paste0('federated/sensitivities/',c('h24','h48','onset','support_only','low_flow','carry0','carry12','raw_support'),'/outcomes/full/model_metadata.csv'))
make_return<-function(site,consensus='frozen_test'){
 p<-file.path(base,site);dir.create(p)
 for(f in required){dir.create(dirname(file.path(p,f)),recursive=TRUE,showWarnings=FALSE);readr::write_csv(tibble::tibble(site=site,status='captured'),file.path(p,f))}
 refresh(p,site,consensus);p
}
refresh<-function(p,site,consensus='frozen_test'){
 f<-list.files(p,recursive=TRUE,full.names=FALSE);f<-setdiff(f,'site_return_manifest.json')
 jsonlite::write_json(list(site=site,primary_consensus_id=consensus,files=as.list(setNames(unname(tools::md5sum(file.path(p,f))),f))),file.path(p,'site_return_manifest.json'),auto_unbox=TRUE)
}
a<-make_return('A');b<-make_return('B')
run<-function(label,paths){system2(file.path(R.home('bin'),'Rscript'),c('diagnostics/site_return/collect.R',shQuote(file.path(base,label)),shQuote(paths)),stdout=file.path(base,paste0(label,'.log')),stderr=file.path(base,paste0(label,'.log')))}
stopifnot(run('valid',c(a,b))==0)
# An empty site CSV can infer character columns where another site has numeric
# values. Empty rows must not force a type conflict or fabricate zero counts.
for(item in list(list(p=a,site='A',n=0L),list(p=b,site='B',n=1L))){
 readr::write_csv(tibble::tibble(site=rep(item$site,item$n),rows=rep(2,item$n)),file.path(item$p,'supplement/optional_empty.csv'))
 refresh(item$p,item$site)
}
stopifnot(run('empty_optional_site',c(a,b))==0)
stopifnot(run('duplicate',c(a,a))!=0)
refresh(b,'B','different_consensus');stopifnot(run('consensus',c(a,b))!=0);refresh(b,'B')
f<-file.path(b,'supplement/sofa_original_components.csv');cat('\n',file=f,append=TRUE)
stopifnot(run('tamper',c(a,b))!=0);refresh(b,'B')
f<-file.path(b,'federated/primary/outcomes/full/model_metadata.csv');readr::write_csv(tibble::tibble(site='B',patient_id='do_not_export'),f);refresh(b,'B')
stopifnot(run('identifiers',c(a,b))!=0)
# New adapted-score phases must be consistent across sites and score specs.
readr::write_csv(tibble::tibble(site='B',status='captured'),f);refresh(b,'B')
setflag<-function(p){m<-jsonlite::fromJSON(file.path(p,'site_return_manifest.json'),simplifyVector=FALSE);m$apache_adapted_extension_complete<-TRUE;jsonlite::write_json(m,file.path(p,'site_return_manifest.json'),auto_unbox=TRUE)}
setflag(a);stopifnot(run('mixed_apache',c(a,b))!=0);refresh(a,'A')
add_apache<-function(p,site,spec='same_score'){
 files<-c('component_completeness.csv','component_counts.csv','cohort_counts.csv','provenance_counts.csv','missing_patterns.csv','source_metadata_quality.csv','model_metadata.csv','transport_metrics.csv','paired_increments.csv','decision_curves.csv','calibration_bins.csv','evaluation_metadata.csv')
 dir<-file.path(p,'federated/apache_adapted');dir.create(dir,recursive=TRUE,showWarnings=FALSE)
 for(k in files)readr::write_csv(tibble::tibble(site=site,status='captured'),file.path(dir,k))
 jsonlite::write_json(list(spec_hash=spec,primary_consensus_id='frozen_test'),file.path(dir,'component_manifest.json'),auto_unbox=TRUE)
 jsonlite::write_json(list(status='passed'),file.path(dir,'audit.json'),auto_unbox=TRUE)
 refresh(p,site);setflag(p)
}
add_apache(a,'A');add_apache(b,'B');stopifnot(run('valid_apache',c(a,b))==0)
add_apache(b,'B','incompatible_score');stopifnot(run('incompatible_apache',c(a,b))!=0)
add_apache(b,'B')
add_curve<-function(p,site){
 f<-file.path(p,'federated/primary/outcomes/full/cumulative_incidence_daily.csv')
 readr::write_csv(tibble::tibble(site=site,consensus_id='frozen_test',phenotype='P1',day=0,event='death',estimate=0,cell_status='observed'),f)
 refresh(p,site);setflag(p)
 m<-jsonlite::fromJSON(file.path(p,'site_return_manifest.json'),simplifyVector=FALSE);m$daily_competing_event_curves_complete<-TRUE
 jsonlite::write_json(m,file.path(p,'site_return_manifest.json'),auto_unbox=TRUE)
}
add_curve(a,'A');stopifnot(run('mixed_curves',c(a,b))!=0)
add_curve(b,'B');stopifnot(run('valid_curves',c(a,b))==0)
add_survival<-function(p,site){
 dir<-file.path(p,'federated/survival_tdauc');dir.create(dir,recursive=TRUE)
 for(f in c('survival_curves.csv','survival_denominators.csv','tdauc.csv','tdroc.csv'))readr::write_csv(tibble::tibble(site=site,status='observed'),file.path(dir,f))
 jsonlite::write_json(list(status='passed',site=site,primary_consensus_id='frozen_test'),file.path(dir,'manifest.json'),auto_unbox=TRUE)
 refresh(p,site);setflag(p)
 m<-jsonlite::fromJSON(file.path(p,'site_return_manifest.json'),simplifyVector=FALSE);m$daily_competing_event_curves_complete<-TRUE;m$survival_tdauc_complete<-TRUE
 jsonlite::write_json(m,file.path(p,'site_return_manifest.json'),auto_unbox=TRUE)
}
add_survival(a,'A');stopifnot(run('mixed_survival',c(a,b))!=0)
add_survival(b,'B');stopifnot(run('valid_survival',c(a,b))==0)
add_pool<-function(p,site,grid='shared_grid'){
 root<-file.path(p,'federated/pooled_survival')
 for(phase in c('times','risksets')){
  d<-file.path(root,phase);dir.create(d,recursive=TRUE,showWarnings=FALSE)
  if(phase=='times'){
   readr::write_csv(tibble::tibble(site=site,endpoint=c('death','death_or_hospice','live_discharge','any_exit')),file.path(d,'cohort_counts.csv'))
   readr::write_csv(tibble::tibble(site=site,events=1L),file.path(d,'event_duration_counts.csv'))
  }else{
   readr::write_csv(tibble::tibble(site=site,n_risk=2L,n_event=1L),file.path(d,'risksets.csv'))
   readr::write_csv(tibble::tibble(site=site,replicate=1L,n_risk=2L,n_event=1L),file.path(d,'bootstrap_risksets.csv'))
  }
  f<-list.files(d,pattern='[.]csv$',full.names=FALSE)
  jsonlite::write_json(list(site=site,status='passed',consensus_id='frozen_test',grid_md5=grid,bootstrap=200,local_AJ_reconstruction_passed=TRUE,local_KM_reconstruction_passed=TRUE,files=as.list(setNames(unname(tools::md5sum(file.path(d,f))),f))),file.path(d,'manifest.json'),auto_unbox=TRUE)
 }
 refresh(p,site)
 m<-jsonlite::fromJSON(file.path(p,'site_return_manifest.json'),simplifyVector=FALSE)
 m$pooled_survival_discharge_complete<-TRUE
 jsonlite::write_json(m,file.path(p,'site_return_manifest.json'),auto_unbox=TRUE)
}
add_pool(a,'A');stopifnot(run('mixed_pool',c(a,b))!=0)
add_pool(b,'B');stopifnot(run('valid_pool',c(a,b))==0)
add_pool(b,'B','different_grid');stopifnot(run('mismatched_pool_grid',c(a,b))!=0)
unlink(base,recursive=TRUE)
cat('PASS: valid collection, duplicate site, consensus mismatch, tamper and nested identifier rejection; adapted-APACHE, daily-curve and survival/tdAUC phase checks\n')
