# Rush-only local audit. The share folder contains aggregate evidence only.
source('code/lib/load.R')
out<-'output/run_RUSH_respiratory_gap_20261030'
if(dir.exists(out))stop('Choose a fresh Rush audit output')
home<-'output/run_federated_rehearsal_20261026/sites/RUSH/local'
a<-readRDS(file.path(home,'discovery.rds'));input<-jsonlite::fromJSON(file.path(a$input_root,'local/input_manifest.json'))
p<-readRDS(file.path(input$baseline_root,'local/prepared.rds'));d<-readRDS(file.path(home,'full/outcomes.rds'))$data
stopifnot(nrow(d)==969L,all(d$time_from_landmark>0))
u<-d[d$hours_UNKNOWN==72,];stopifnot(nrow(u)==126L,all(u$assignment_status=='insufficient_support'))
keys<-d |> dplyr::select(hospitalization_id,t0,admission_year,hours_UNKNOWN)
window<-function(x,ts='recorded_dttm')x |> dplyr::inner_join(keys,by='hospitalization_id') |> dplyr::filter(.data[[ts]]>=t0,.data[[ts]]<t0+72*3600)
r<-p$tables$respiratory_support;rw<-window(r)
uw<-rw[rw$hospitalization_id%in%u$hospitalization_id,]
ids_any<-unique(uw$hospitalization_id)
vu<-window(p$tables$vitals) |> dplyr::filter(hospitalization_id%in%u$hospitalization_id,is.finite(num(vital_value)))
lu<-window(p$tables$labs,'lab_result_dttm') |> dplyr::filter(hospitalization_id%in%u$hospitalization_id,is.finite(num(lab_value_numeric)))
sum_ids<-function(x)length(intersect(u$hospitalization_id,unique(x)))
summary<-tibble::tibble(metric=c('landmark_hospitalizations','any_known_support_hours','all72h_unknown_support','no_resp_rows_first72h_among_all_unknown','any_resp_rows_first72h_among_all_unknown','no_resp_rows_entire_hospitalization_among_all_unknown','any_fio2_first72h_among_all_unknown','positive_flow_first72h_among_all_unknown','any_numeric_vital_first72h_among_all_unknown','any_spo2_first72h_among_all_unknown','any_numeric_lab_first72h_among_all_unknown','unknown_support_hours','all_support_hours'),
 n=c(nrow(d),sum(d$hours_UNKNOWN<72),nrow(u),nrow(u)-length(ids_any),length(ids_any),nrow(u)-sum_ids(r$hospitalization_id),sum_ids(uw$hospitalization_id[which(!is.na(num(uw$fio2_set)))]),sum_ids(uw$hospitalization_id[which(num(uw$lpm_set)>0)]),sum_ids(vu$hospitalization_id),sum_ids(vu$hospitalization_id[tolower(vu$vital_category)=='spo2']),sum_ids(lu$hospitalization_id),sum(d$hours_UNKNOWN),72*nrow(d)))
year<-d |> dplyr::group_by(admission_year) |> dplyr::summarise(n=dplyr::n(),all72h_unknown_n=sum(hours_UNKNOWN==72),known_some_hours_n=sum(hours_UNKNOWN<72),unknown_support_hours=sum(hours_UNKNOWN),all_support_hours=72*dplyr::n(),.groups='drop')
year$all72h_unknown_pct<-100*year$all72h_unknown_n/year$n
# ICU exposure: union of closed, documented ICU ADT intervals only; open-ended
# intervals are flagged separately rather than assumed to cover the whole window.
ad<-p$tables$adt |> dplyr::inner_join(u |> dplyr::select(hospitalization_id,t0),by='hospitalization_id') |> dplyr::filter(tolower(location_category)=='icu',in_dttm<t0+72*3600,is.na(out_dttm)|out_dttm>t0)
exposure<-vapply(u$hospitalization_id,function(id){
 z<-ad[ad$hospitalization_id==id&!is.na(ad$out_dttm),]
 if(!nrow(z))return(0)
 lo<-pmax(0,hours_between(z$in_dttm,z$t0));hi<-pmin(72,hours_between(z$out_dttm,z$t0));ord<-order(lo);lo<-lo[ord];hi<-hi[ord]
 total<-0;end<-0
 for(j in seq_along(lo))if(hi[j]>lo[j]){total<-total+max(0,hi[j]-max(end,lo[j]));end<-max(end,hi[j])}
 total
},numeric(1))
icu<-tibble::tibble(metric=c('unknown_group_n','closed_ICU_interval_hours_median','closed_ICU_interval_hours_q1','closed_ICU_interval_hours_q3','at_least24_closed_ICU_hours_n','at_least72_closed_ICU_hours_n','any_open_ended_ICU_interval_n'),
 value=c(nrow(u),median(exposure),unname(quantile(exposure,.25)),unname(quantile(exposure,.75)),sum(exposure>=24),sum(exposure>=72),dplyr::n_distinct(ad$hospitalization_id[is.na(ad$out_dttm)])))
source_labels<-uw |> dplyr::count(device_category,name='rows')
dir.create(file.path(out,'share'),recursive=TRUE)
write<-function(x,name){stopifnot(!any(names(x)%in%c('patient_id','hospitalization_id','t0','recorded_dttm','in_dttm','out_dttm')));readr::write_csv(x,file.path(out,'share',name),na='')}
write(summary,'respiratory_gap_summary.csv');write(year,'gap_by_admission_year.csv');write(icu,'icu_exposure_summary.csv');write(source_labels,'source_device_labels_in_unknown_group.csv')
stopifnot(summary$n[summary$metric=='no_resp_rows_first72h_among_all_unknown']==125L,sum(year$n)==969L,sum(year$all72h_unknown_n)==126L)
files<-list.files(file.path(out,'share'),full.names=TRUE)
jsonlite::write_json(list(status='passed',site='RUSH',scope='72h landmark cohort, frozen three-site primary analysis; no support-coverage exclusion',
 respiratory_processing='CLIF waterfall support, source-evidence ARF; first72h after first ICU arrival',
 source_release='Analyzed Rush CLIF2.1 extract cached in October22 run; not a fresh audit of a newer release',
 conclusion='Respiratory-table coverage gap; extraction versus source-documentation cause unconfirmed',
 source_md5=as.list(tools::md5sum(c(file.path(input$baseline_root,'local/prepared.rds'),file.path(home,'full/outcomes.rds'),'diagnostics/rush_respiratory_gap/audit.R'))),
 files=as.list(setNames(unname(tools::md5sum(files)),basename(files)))),file.path(out,'audit_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
print(summary);print(year);print(icu);print(source_labels)
