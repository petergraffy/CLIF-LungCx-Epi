# Exact aggregate extension of frozen clinical/model caches. No model refits.
source('diagnostics/apache_adapted/lib.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=2)stop('Usage: site_counts.R SITE NEW_OUTPUT')
site<-a[1];out<-a[2];if(dir.exists(out))stop('Choose fresh exact-count output')
stopifnot(aggregate_exact());dir.create(out,recursive=TRUE)
primary<-ap_primary();revision<-Sys.getenv('CLIF_REVISION_OUTPUT',unset='output/run_revision_gaps_20261026')
home<-file.path(primary,'sites',site);disc<-readRDS(file.path(home,'local/discovery.rds'))
input<-jsonlite::fromJSON(file.path(disc$input_root,'local/input_manifest.json'));p<-readRDS(file.path(input$baseline_root,'local/prepared.rds'))
cfg<-list(site_name=site,analysis=p$parameters);cfg$analysis$min_cell<-0
write<-function(x,name){
 x<-promote_exact_counts(x)
 stopifnot(!any(names(x)%in%c('patient_id','hospitalization_id','t0','event_time','recorded_dttm','admission_dttm','discharge_dttm','sequence','marker','threshold')))
 if(!'site'%in%names(x))x$site<-rep(site,nrow(x))
 dest<-file.path(out,paste0(name,'.csv'));dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE);readr::write_csv(x,dest,na='')
}
flow<-landmark_cohort(p$cohort,72) |> dplyr::count(landmark_reason,name='n');write(flow,'cohort/landmark_flow')
write(p$cohort |> dplyr::count(admission_year,event_type,valid_outcome,name='n'),'cohort/annual_outcomes')
rebuilt<-build_cohort(p$tables,cfg);stopifnot(setequal(rebuilt$cohort$hospitalization_id,p$cohort$hospitalization_id));write(rebuilt$flow,'cohort/screening_flow')
counts<-list();models<-list()
capture<-function(d,analysis,version,fit=NULL){
 groups<-if('phenotype'%in%names(d))c('Overall',sort(unique(d$phenotype))) else 'Overall'
 for(g in groups){z<-if(g=='Overall')d else d[d$phenotype==g,]
  t<-if(analysis!='absorbing'&&'time_from_landmark'%in%names(z))z$time_from_landmark else hours_between(z$event_time,z$t0)/24
  counts[[paste(analysis,version,g)]]<<-tibble::tibble(analysis=analysis,version=version,group=g,n=nrow(z),site_distinct_patients=dplyr::n_distinct(z$patient_id),unknown_outcomes=sum(!z$valid_outcome),death28=sum(z$valid_outcome&t<=28&z$event_type=='death'),hospice28=sum(z$valid_outcome&t<=28&z$event_type=='hospice'),live_discharge28=sum(z$valid_outcome&t<=28&z$event_type=='live_discharge'),no_terminal_event28=sum(z$valid_outcome&t>28))
 }
 if(!is.null(fit))for(k in names(fit)){
  f<-fit[[k]];ok<-inherits(f,'coxph')
  models[[paste(analysis,version,k)]]<<-tibble::tibble(analysis=analysis,version=version,model=k,n=if(ok)f$n else NA_integer_,events=if(ok)f$nevent else NA_integer_,status=if(ok)'saved_model' else 'not_estimable')
 }
}
capture(p$cohort,'eligible','t0')
for(v in c('full',paste0('loo_',site))){x<-readRDS(file.path(home,'local',v,'outcomes.rds'));capture(x$data,'primary',v,x$cox)}
for(variant in c('h24','h48','onset','support_only','low_flow','carry0','carry12','raw_support'))for(v in c('full',paste0('loo_',site))){
 x<-readRDS(file.path(revision,'sensitivities',variant,'sites',site,'local',v,'outcomes.rds'));capture(x$data,variant,v,x$cox)
}
x<-readRDS(file.path(revision,'sensitivities/absorbing_v2/sites',site,'local/absorbing_summary.rds'))
capture(x$data,'absorbing','descriptive')
write(dplyr::bind_rows(counts),'cohort/all_analysis_counts');write(dplyr::bind_rows(models),'models/cox_denominators')
for(pair in list(c(file.path(revision,'sites',site,'local/model_denominators.csv'),'models/reviewer_denominators'),c(file.path(ap_root(),'sites',site,'local/model_denominators_exact.csv'),'models/apache_training_denominators'))){write(readr::read_csv(pair[1],show_col_types=FALSE),pair[2])}
ap<-readRDS(file.path(ap_root(),'sites',site,'local/evaluation.rds'));ev<-list()
for(k in names(ap)){x<-ap[[k]];ev[[k]]<-tibble::tibble(cohort_endpoint=k,n=x$n,events=x$events,non_events=x$n-x$events,patients=dplyr::n_distinct(x$patient_id))}
write(dplyr::bind_rows(ev),'models/apache_evaluation_denominators')
so<-readRDS(file.path(revision,'sites',site,'local/score_transport.rds'));ev<-list()
for(k in names(so)){x<-so[[k]];ev[[k]]<-tibble::tibble(cohort_endpoint=k,n=x$n,events=x$events,non_events=x$n-x$events,patients=dplyr::n_distinct(x$patient_id))}
write(dplyr::bind_rows(ev),'models/sofa_evaluation_denominators')
apc<-readRDS(file.path(ap_root(),'sites',site,'local/components.rds'))
cc<-lapply(c(APACHE_COMPONENTS,'adapted_total','observed_hct_total'),function(k){v<-apc[[paste0('apache_',k)]];tibble::tibble(component=k,n=nrow(apc),observed_n=sum(is.finite(v)),missing_n=sum(!is.finite(v)))})
write(dplyr::bind_rows(cc),'severity/apache_component_counts')
patterns<-tibble::tibble(pattern=apply(is.na(apc[paste0('apache_',APACHE_COMPONENTS)]),1,paste0,collapse='')) |> dplyr::count(pattern,name='n');write(patterns,'severity/apache_missing_patterns')
full<-readRDS(file.path(home,'local/full/outcomes.rds'))$data
panel<-disc$panel |> dplyr::left_join(full |> dplyr::select(hospitalization_id,phenotype),by='hospitalization_id')
for(domain in c('support','hypox','hyper','joint')){
 z<-panel |> dplyr::count(phenotype,h,state=.data[[domain]],name='n') |> dplyr::group_by(phenotype,h) |> dplyr::mutate(denominator=sum(n),proportion=n/denominator) |> dplyr::ungroup()
 write(z,paste0('profiles/hourly_',domain))
}
write(disc$data |> dplyr::count(local_cluster,name='n') |> dplyr::left_join(disc$stability,by='local_cluster'),'discovery/local_cluster_counts')
# These are already aggregate diagnostics, saved before historical masking.
exact<-file.path('output/run_site_returns_20261026',site,'local/supplement_exact')
if(dir.exists(exact))for(f in list.files(exact,pattern='[.]csv$',full.names=TRUE)){
 if(basename(f)%in%c('return_file_inventory.csv','analysis_return_inventory.csv'))next
 write(readr::read_csv(f,show_col_types=FALSE),paste0('supplement/',tools::file_path_sans_ext(basename(f))))
}
files<-list.files(out,recursive=TRUE,full.names=FALSE)
jsonlite::write_json(list(status='passed',site=site,aggregate_policy='exact counts; no disclosure suppression',primary_consensus_id=jsonlite::fromJSON(file.path(primary,'coordinator/full/frozen_manifest.json'))$consensus_id,
 original_models_refitted=FALSE,patient_level_identifiers_returned=FALSE,patient_counts_not_globally_additive=TRUE,files=as.list(setNames(unname(tools::md5sum(file.path(out,files))),files))),file.path(out,'manifest.json'),auto_unbox=TRUE,pretty=TRUE)
cat(site,'exact cohort/model/profile/supplement counts exported\n')
