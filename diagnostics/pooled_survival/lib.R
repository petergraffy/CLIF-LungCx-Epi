source('code/lib/load.R')
PS_TICKS_PER_DAY<-86400*1e6
ps_tick<-function(time)round(time*PS_TICKS_PER_DAY)
ps_time<-function(tick)tick/PS_TICKS_PER_DAY
ps_data<-function(site,primary){
 home<-file.path(primary,'sites',site);dis<-readRDS(file.path(home,'local/discovery.rds'))
 im<-jsonlite::fromJSON(file.path(dis$input_root,'local/input_manifest.json'))
 stopifnot(im$site==site)
 p<-readRDS(file.path(im$baseline_root,'local/prepared.rds'))
 full<-readRDS(file.path(home,'local/full/outcomes.rds'))$data
 all<-p$cohort;all$time<-hours_between(all$event_time,all$t0)/24
 full$time<-full$time_from_landmark
 all$time<-ps_time(ps_tick(all$time));full$time<-ps_time(ps_tick(full$time))
 stopifnot(all(is.finite(all$time)&all$time>0),all(is.finite(full$time)&full$time>0))
 sets<-list(icu_admission_Overall=all,landmark72h_Overall=full)
 phenotypes<-rownames(readRDS(file.path(primary,'coordinator/full/prototypes.rds')))
 for(ph in c(phenotypes,'I'))sets[[paste0('landmark72h_',ph)]]<-full[full$phenotype==ph,]
 list(sets=sets,patients=unique(all$patient_id))
}
ps_event<-function(d,endpoint){
 types<-switch(endpoint,death='death',death_or_hospice=c('death','hospice'),live_discharge='live_discharge',any_exit=c('death','hospice','live_discharge'),stop('Unknown endpoint'))
 as.integer(d$valid_outcome & d$event_type%in%types)
}
# Aalen-Johansen incidence: competing exits remove probability mass permanently.
ps_incidence<-function(risk,events,all_events){
 stopifnot(length(risk)==length(events),length(events)==length(all_events),all(events>=0),all(events<=all_events),all(all_events<=risk))
 survival<-ps_product(risk,all_events)
 out<-cumsum(c(1,head(survival,-1))*ifelse(risk>0,events/pmax(risk,1),0))
 if(risk[1]==0)out[]<-NA_real_
 out
}
ps_export<-function(x,p){
 stopifnot(!any(names(x)%in%c('patient_id','hospitalization_id','t0','event_time','recorded_dttm','admission_dttm','discharge_dttm','sequence')))
 dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);readr::write_csv(x,p,na='')
}
ps_counts<-function(T,event,grid,weights=rep(1,length(T))){
 # Event/censor ties: subjects remain at risk for events at their exit time.
 ticks<-ps_tick(T);g_ticks<-ps_tick(grid)
 order<-order(ticks);cum<-c(0,cumsum(weights[order]));less<-findInterval(g_ticks,ticks[order],left.open=TRUE)
 risk<-sum(weights)-cum[less+1L]
 ev<-which(event==1L&T<=28);index<-match(ticks[ev],g_ticks);stopifnot(!anyNA(index),!anyDuplicated(g_ticks))
 events<-numeric(length(grid));if(length(ev)){
  q<-rowsum(weights[ev],index,reorder=FALSE);events[as.integer(rownames(q))]<-q[,1]
 }
 stopifnot(all(events<=risk),all(risk>=0),all(events>=0))
 list(risk=as.numeric(risk),events=events)
}
ps_product<-function(risk,events){
 stopifnot(all(risk>=0),all(events>=0),all(events<=risk))
 factor<-ifelse(risk>0,1-events/pmax(risk,1),1)
 value<-cumprod(factor)
 if(risk[1]==0)value[]<-NA_real_
 value
}
