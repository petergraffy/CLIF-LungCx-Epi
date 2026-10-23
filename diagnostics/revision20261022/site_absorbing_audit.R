source('diagnostics/revision20261022/absorbing_contract.R')
site<-commandArgs(trailingOnly=TRUE)[1];home<-file.path(fed_root(),'sites',site);a<-readRDS(file.path(home,'local/discovery.rds'))
d<-a$data;p<-a$panel |> dplyr::left_join(d |> dplyr::select(hospitalization_id,t0,event_time,event_type),by='hospitalization_id')
expected<-floor(hours_between(p$event_time,p$t0))<=p$h
stopifnot(all(p$terminal==expected));names<-c(death='DEATH',hospice='HOSPICE',live_discharge='LIVE_DISCHARGE');stopifnot(all(p$support[expected]==names[p$event_type[expected]]))
mo<-fed_modal(a$panel);stopifnot(!anyNA(mo$state),setequal(unique(mo$bin[mo$domain=='support']),0:11))
z<-readRDS(file.path(home,'local/absorbing_summary.rds'));stopifnot(nrow(z$data)==nrow(d),!anyNA(z$data$descriptive_group),!dir.exists(file.path(home,'export/outcomes')))
fed_json(list(status='passed',site=site,n=nrow(d),terminal_states_start_at_event_hour=TRUE,terminal_transitions_not_clinically_ranked=TRUE,no_outcome_regression=TRUE,all_cases_retained=TRUE),file.path(home,'export/absorbing_audit.json'))
