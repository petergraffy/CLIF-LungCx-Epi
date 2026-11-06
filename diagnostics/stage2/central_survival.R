# Convert exact aggregate exit histograms to the existing pooled AJ/KM interface.
source('scripts/deployment_lib.R');source('diagnostics/stage2/lib.R');source('diagnostics/pooled_survival/lib.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)<2L)stop('Usage: central_survival.R NEW_ROOT SITE_SURVIVAL [...]')
out<-a[1];shares<-a[-1];if(dir.exists(out))stop('Choose fresh central survival root')
metas<-lapply(shares,deploy_validate_inventory,'manifest.json');sites<-vapply(metas,`[[`,character(1),'site');B<-metas[[1]]$bootstrap
stopifnot(!anyDuplicated(sites),length(unique(vapply(metas,`[[`,character(1),'consensus_id')))==1L,all(vapply(metas,function(m)m$bootstrap==B,logical(1))))
hist<-lapply(shares,function(p)readr::read_csv(file.path(p,'exit_histograms.csv'),show_col_types=FALSE))
counts<-lapply(shares,function(p)readr::read_csv(file.path(p,'cohort_counts.csv'),show_col_types=FALSE));curves<-unique(counts[[1]]$curve)
stopifnot(all(vapply(counts,function(z)setequal(z$curve,curves),logical(1))))
grids<-list()
for(curve in curves){ticks<-sort(unique(c(ps_tick(c(0,7,14,21,28)),unlist(lapply(hist,function(z)z$time_tick[z$curve==curve])))))
 for(endpoint in c('death','death_or_hospice','live_discharge','any_exit'))grids[[paste(curve,endpoint)]]<-tibble::tibble(curve=curve,endpoint=endpoint,time_tick=ticks,time=ps_time(ticks))
}
grid<-dplyr::bind_rows(grids);grid$grid_id<-seq_len(nrow(grid));center<-file.path(out,'grid');ps_export(grid,file.path(center,'shared_grid.csv'))
ghash<-unname(tools::md5sum(file.path(center,'shared_grid.csv')))
deploy_json(list(status='passed',sites=sites,consensus_id=metas[[1]]$consensus_id,grid_md5=ghash,patient_level_inputs_read=FALSE),file.path(center,'manifest.json'))
for(i in seq_along(sites)){
 nominal<-list();boots<-list()
 for(curve in curves){
  z<-hist[[i]][hist[[i]]$curve==curve,];ct<-counts[[i]][counts[[i]]$curve==curve,];stopifnot(nrow(ct)==1L,sum(z$n)==ct$n,all(z$n>=0),all(z$site==sites[i]))
  for(endpoint in c('death','death_or_hospice','live_discharge','any_exit')){
   g<-grid[grid$curve==curve & grid$endpoint==endpoint,];types<-switch(endpoint,death='death',death_or_hospice=c('death','hospice'),live_discharge='live_discharge',any_exit=c('death','hospice','live_discharge'))
   isevent<-z$exit_type%in%types;idx<-match(z$time_tick,g$time_tick);stopifnot(!anyNA(idx))
   count_events<-function(w){r<-numeric(nrow(g));if(any(isevent)){q<-rowsum(w[isevent],idx[isevent],reorder=FALSE);r[as.integer(rownames(q))]<-q[,1]};r}
   risk<-s2_hist_risk(z$time_tick,z$n,g$time_tick,ct$n)
   nominal[[paste(curve,endpoint)]]<-g |> dplyr::mutate(site=sites[i],n_risk=risk,n_event=count_events(z$n),n=ct$n,max_followup=ct$max_followup)
   for(b in seq_len(B)){
    w<-z[[sprintf('n_b%03d',b)]];stopifnot(all(w>=0))
    boots[[paste(curve,endpoint,b)]]<-tibble::tibble(site=sites[i],grid_id=g$grid_id,replicate=b,n_risk=s2_hist_risk(z$time_tick,w,g$time_tick,sum(w)),n_event=count_events(w),max_followup=if(any(w>0))ps_time(max(z$time_tick[w>0])) else 0)
   }
  }
 }
 dir<-file.path(out,'sites',sites[i]);ps_export(dplyr::bind_rows(nominal),file.path(dir,'risksets.csv'));ps_export(dplyr::bind_rows(boots),file.path(dir,'bootstrap_risksets.csv'));f<-list.files(dir,full.names=TRUE)
 deploy_json(list(status='passed',site=sites[i],consensus_id=metas[[i]]$consensus_id,bootstrap=B,grid_md5=ghash,local_KM_reconstruction_passed=TRUE,local_AJ_reconstruction_passed=TRUE,reconstruction='from exact aggregate exit histograms; all-time risk-set identity checked by one-pass regression tests',files=as.list(setNames(unname(tools::md5sum(f)),basename(f)))),file.path(dir,'manifest.json'))
}
