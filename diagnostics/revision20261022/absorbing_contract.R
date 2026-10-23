source('diagnostics/revision20261022/lib.R')
REV_WINDOW<-72L
base_spec<-fed_spec()
fed_spec<-function()modifyList(base_spec,list(schema='revision-absorbing-v1',exposure='waterfall support with raw ARF and absorbing death/hospice/live-discharge; descriptive only'))
fed_domains<-function()list(joint=c(sort(as.vector(outer(SUPPORT_STATES,PHYS_STATES,paste,sep='|'))),TERMINAL_STATES),support=c(SUPPORT_STATES,TERMINAL_STATES),hypox=c('H0','H1','HU'),hyper=c('C0','C1','CU'),change=c('down','same','up','unknown'))
# Changes into a terminal event are descriptive, never an escalation rank.
original_modal<-fed_modal
fed_modal<-function(panel){
 p<-original_modal(panel)
 q<-p |> dplyr::filter(domain=='support') |> dplyr::arrange(hospitalization_id,bin) |> dplyr::group_by(hospitalization_id) |> dplyr::mutate(terminal_transition=state%in%TERMINAL_STATES|dplyr::lag(state)%in%TERMINAL_STATES) |> dplyr::ungroup()
 terminal_keys<-paste(q$hospitalization_id[q$terminal_transition],q$bin[q$terminal_transition])
 ix<-p$domain=='change' & (paste(p$hospitalization_id,p$bin)%in%terminal_keys|is.na(p$state));p$state[ix]<-'unknown';p
}
revision_input<-function(site){
 a<-readRDS(paste0('output/run_',site,'_waterfall_20261022/local/waterfall_derivation.rds'))$analyses$h72_absorbing
 ix<-!a$panel$terminal;a$panel$physiology[ix]<-a$panel$old_physiology[ix];a$panel$hypox[ix]<-sub('_.*','',a$panel$physiology[ix]);a$panel$hyper[ix]<-sub('.*_','',a$panel$physiology[ix]);a$panel$joint[ix]<-paste(a$panel$support[ix],a$panel$physiology[ix],sep='|')
 a$data<-a$data |> dplyr::select(-dplyr::any_of(c('phenotype','local_cluster')));a
}
