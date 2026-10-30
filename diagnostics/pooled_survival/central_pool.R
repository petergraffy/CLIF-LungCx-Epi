# All inputs here are aggregate site returns; no clinical cache access.
source('diagnostics/pooled_survival/lib.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)<3)stop('Usage: central_pool.R NEW_OUTPUT CENTRAL_GRID SITE_RISK_SHARE [...]')
out<-a[1];center<-a[2];shares<-a[-c(1,2)];if(dir.exists(out))stop('Choose fresh pooled output')
cm<-jsonlite::fromJSON(file.path(center,'manifest.json'));grid<-readr::read_csv(file.path(center,'shared_grid.csv'),show_col_types=FALSE)
stopifnot(unname(tools::md5sum(file.path(center,'shared_grid.csv')))==cm$grid_md5)
metas<-lapply(shares,function(p)jsonlite::fromJSON(file.path(p,'manifest.json')))
sites<-vapply(metas,function(m)m$site,character(1));B<-200L;L<-nrow(grid)
stopifnot(!anyDuplicated(sites),setequal(sites,cm$sites),all(vapply(metas,function(m)m$grid_md5==cm$grid_md5&&m$consensus_id==cm$consensus_id&&m$bootstrap==B&&m$local_KM_reconstruction_passed&&m$local_AJ_reconstruction_passed,logical(1))))
nominal<-list();risk<-events<-maxfu<-matrix(0,L,B)
for(i in seq_along(shares)){
 for(k in names(metas[[i]]$files))stopifnot(unname(tools::md5sum(file.path(shares[i],k)))==metas[[i]]$files[[k]])
 z<-readr::read_csv(file.path(shares[i],'risksets.csv'),show_col_types=FALSE);stopifnot(nrow(z)==L,all(z$site==sites[i]),setequal(z$grid_id,grid$grid_id))
 nominal[[i]]<-z
 b<-readr::read_csv(file.path(shares[i],'bootstrap_risksets.csv'),show_col_types=FALSE)
 index<-b$grid_id+(b$replicate-1L)*L
 stopifnot(nrow(b)==L*B,!anyDuplicated(index),all(b$site==sites[i]),all(index>=1&index<=L*B))
 risk[index]<-risk[index]+b$n_risk;events[index]<-events[index]+b$n_event;maxfu[index]<-pmax(maxfu[index],b$max_followup)
 rm(b,index);invisible(gc())
}
site<-dplyr::bind_rows(nominal)
nom<-site |> dplyr::group_by(grid_id,curve,endpoint,time) |> dplyr::summarise(n_risk=sum(n_risk),n_event=sum(n_event),contributing_sites=sum(n>0),n=sum(n),max_followup=max(max_followup),.groups='drop') |> dplyr::arrange(grid_id)
pooled<-list();site_curves<-list();audit<-list()
keys<-unique(grid[c('curve','endpoint')])
for(i in seq_len(nrow(keys))){
 k<-keys[i,];z<-nom[nom$curve==k$curve&nom$endpoint==k$endpoint,] |> dplyr::arrange(time)
 ids<-z$grid_id;point<-ps_product(z$n_risk,z$n_event);point[z$time>z$max_followup]<-NA_real_
 bs<-vapply(seq_len(B),function(b){v<-ps_product(risk[ids,b],events[ids,b]);v[z$time>maxfu[ids,b]]<-NA_real_;v},numeric(length(ids)))
 finite<-rowSums(is.finite(bs))
 ci<-t(vapply(seq_len(nrow(bs)),function(j)if(finite[j]>=100)unname(quantile(bs[j,],c(.025,.975),na.rm=TRUE)) else c(NA_real_,NA_real_),numeric(2)))
 z$survival<-point;z$conf_low<-ci[,1];z$conf_high<-ci[,2];z$bootstrap_replicates<-finite
 z$clock<-if(grepl('^icu_admission_',k$curve))'icu_admission' else'landmark72h';z$group<-sub('^(icu_admission_|landmark72h_)','',k$curve)
 z$status<-ifelse(z$time>z$max_followup,'beyond_observed_followup',ifelse(is.finite(point),'observed','no_hospitalizations'))
 pooled[[i]]<-z
 stopifnot(all(diff(point[is.finite(point)])<=1e-12),z$n_risk[1]==z$n[1])
 if(z$n[1]>0)stopifnot(point[1]==1) else stopifnot(all(is.na(point)))
 audit[[i]]<-tibble::tibble(curve=k$curve,endpoint=k$endpoint,n=z$n[1],at_risk_at_zero=z$n_risk[1],events28=sum(z$n_event),bootstrap=200L)
 for(s in sites){
  q<-site[site$site==s&site$curve==k$curve&site$endpoint==k$endpoint,] |> dplyr::arrange(time)
  q$survival<-ps_product(q$n_risk,q$n_event);q$survival[q$time>q$max_followup]<-NA_real_;site_curves[[paste(i,s)]]<-q
 }
}
# Compute live-discharge incidence from summed counts, never from averaged site curves.
incidence<-list();site_incidence<-list()
for(curve_key in unique(grid$curve)){
 z<-nom[nom$curve==curve_key&nom$endpoint=='live_discharge',] |> dplyr::arrange(time)
 exits<-nom[nom$curve==curve_key&nom$endpoint=='any_exit',] |> dplyr::arrange(time)
 stopifnot(identical(z$time,exits$time),identical(z$n_risk,exits$n_risk))
 point<-ps_incidence(z$n_risk,z$n_event,exits$n_event)
 bs<-vapply(seq_len(B),function(b){
  v<-ps_incidence(risk[z$grid_id,b],events[z$grid_id,b],events[exits$grid_id,b])
  v[z$time>maxfu[z$grid_id,b]]<-NA_real_;v
 },numeric(nrow(z)))
 finite<-rowSums(is.finite(bs))
 ci<-t(vapply(seq_len(nrow(bs)),function(j)if(finite[j]>=100)unname(quantile(bs[j,],c(.025,.975),na.rm=TRUE)) else c(NA_real_,NA_real_),numeric(2)))
 z$estimate<-point;z$estimate[z$time>z$max_followup]<-NA_real_
 z$conf_low<-ci[,1];z$conf_high<-ci[,2];z$bootstrap_replicates<-finite
 z$n_competing<-exits$n_event-z$n_event
 stopifnot(all(diff(point[is.finite(point)])>=-1e-12),all(point[is.finite(point)]<=1+1e-12))
 incidence[[curve_key]]<-z
 for(s in sites){
  q<-site[site$site==s&site$curve==curve_key&site$endpoint=='live_discharge',] |> dplyr::arrange(time)
  e<-site[site$site==s&site$curve==curve_key&site$endpoint=='any_exit',] |> dplyr::arrange(time)
  q$estimate<-ps_incidence(q$n_risk,q$n_event,e$n_event)
  q$estimate[q$time>q$max_followup]<-NA_real_;site_incidence[[paste(curve_key,s)]]<-q
 }
}
ps_export(dplyr::bind_rows(incidence),file.path(out,'pooled_discharge_incidence.csv'))
ps_export(dplyr::bind_rows(site_incidence),file.path(out,'site_discharge_incidence.csv'))
pooled<-dplyr::bind_rows(pooled)
ps_export(pooled,file.path(out,'pooled_survival.csv'))
ps_export(pooled |> dplyr::filter(time%in%c(0,7,14,21,28)),file.path(out,'risk_tables.csv'))
ps_export(dplyr::bind_rows(site_curves),file.path(out,'site_survival.csv'))
ps_export(dplyr::bind_rows(audit),file.path(out,'pooling_denominators.csv'))
jsonlite::write_json(list(status='passed',sites=sites,consensus_id=cm$consensus_id,patient_level_inputs_read=FALSE,method='Pooled product-limit: product over exact shared event-duration grid of1 - sum_site(events)/sum_site(risk)',
 discharge_endpoint='Live non-hospice discharge; death and hospice are competing events; unknown exits censored',discharge_method='Aalen-Johansen cumulative incidence from exact summed risk and all-exit counts; 200 within-site patient-cluster bootstrap intervals',target='Observed combined hospitalization cohort; risk-set weighting, not equal-site averaging',bootstrap='200 independent patient-cluster resamples within site; all admissions carried together; pointwise percentile intervals; conditional on these three sites',
 clocks='All eligible admissions at ICU t0; frozen primary phenotypes at72h',limitations=c('Unadjusted curves; retain site-specific heterogeneity and adjusted meta-analysis','Competing hospital exits censored in cause-specific KM; not post-discharge all-cause survival','Unknown whole-cohort dispositions censored at observed discharge','Cannot link repeated patients across institutions'),
 shared_data='Relative follow-up duration grid and exact aggregate risk/event counts; no identifiers or calendar timestamps',clinical_models_refitted=FALSE,source_md5=as.list(tools::md5sum(c('diagnostics/pooled_survival/lib.R','diagnostics/pooled_survival/central_pool.R')))),file.path(out,'pooling_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
cat('Exact pooled KM curves and patient-cluster bootstrap intervals passed\n')
