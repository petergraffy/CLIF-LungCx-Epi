source('diagnostics/federated/model_lib.R')
root<-Sys.getenv('CLIF_MIMIC_OUTPUT',unset='output/run_MIMIC_external_20261023');local<-file.path(root,'local');dir.create(file.path(root,'aggregate'),recursive=TRUE,showWarnings=FALSE)
path<-file.path(local,'source_tables');cfg<-list(site_name='MIMIC',repo=getwd(),tables_path=path,file_type='parquet',analysis=modifyList(run_defaults(),list(run_id='external_20261023',start_date='1900-01-01',end_date='2300-01-01')))
# Deidentified dates preserve intervals, not the consortium's calendar window.
read_table<-function(nm,cols=NULL,ids=NULL,patients=NULL){
 x<-arrow::open_dataset(file.path(path,paste0('clif_',nm,'.parquet')))
 if(!is.null(cols))x<-dplyr::select(x,dplyr::all_of(intersect(cols,names(x))))
 if(!is.null(ids))x<-dplyr::filter(x,hospitalization_id%in%ids)
 if(!is.null(patients))x<-dplyr::filter(x,patient_id%in%patients)
 z<-dplyr::collect(x)
 for(ts in grep('_dttm$',names(z),value=TRUE))z[[ts]]<-utc(z[[ts]])
 z
}
p<-read_table('patient',c('patient_id','sex_category','race_category','ethnicity_category','death_dttm'))
h<-read_table('hospitalization',c('patient_id','hospitalization_id','admission_dttm','discharge_dttm','age_at_admission','discharge_category'))
adt<-read_table('adt',c('hospitalization_id','in_dttm','out_dttm','location_category'))
dx<-read_table('hospital_diagnosis',c('hospitalization_id','diagnosis_code','diagnosis_code_format','poa_present'))
code<-toupper(gsub('[^A-Za-z0-9]','',dx$diagnosis_code));fmt<-toupper(dx$diagnosis_code_format)
is_lung<-(!is.na(fmt)&grepl('10',fmt)&grepl('^C34',code))|(!is.na(fmt)&grepl('9',fmt)&grepl('^162',code));is_lung[is.na(is_lung)]<-FALSE
lung_ids<-unique(dx$hospitalization_id[is_lung]);poa_known<-sum(!is.na(dx$poa_present))
# This copy is an explicit cohort-selection adapter, NEVER observed POA evidence.
# Raw diagnosis table remains unchanged and is not used for baseline covariates.
eligibility_dx<-dx[is_lung,];eligibility_dx$poa_present<-1L
built<-build_cohort(list(patient=p,hospitalization=h,adt=adt,hospital_diagnosis=eligibility_dx),cfg)
d<-built$cohort;stopifnot(all(d$hospitalization_id%in%lung_ids))
prior<-h[h$hospitalization_id%in%lung_ids & !is.na(h$discharge_dttm),] |> dplyr::group_by(patient_id) |> dplyr::summarise(first_completed_lung_coded_encounter=min(discharge_dttm),.groups='drop')
d<-d |> dplyr::left_join(prior,by='patient_id');d$prior_completed_lung_cancer_encounter<-!is.na(d$first_completed_lung_coded_encounter)&d$first_completed_lung_coded_encounter<=d$admission_dttm
# Do not fabricate POA comorbidity or cancer-stage covariates from discharge codes.
d$charlson_noncancer<-NA_real_;d$advanced_cancer<-NA_real_;d$covid<-NA
ids<-d$hospitalization_id;patients<-unique(d$patient_id)
tables<-list(patient=p[p$patient_id%in%patients,],hospitalization=h[h$hospitalization_id%in%ids,],adt=adt[adt$hospitalization_id%in%ids,],hospital_diagnosis=dx[dx$hospitalization_id%in%ids,])
columns<-list(respiratory_support=NULL,vitals=c('hospitalization_id','recorded_dttm','vital_category','vital_value'),labs=c('hospitalization_id','lab_result_dttm','lab_category','lab_value_numeric'),medication_admin_continuous=c('hospitalization_id','admin_dttm','med_category','med_dose','med_dose_unit'),patient_assessments=c('hospitalization_id','recorded_dttm','assessment_category','numerical_value'))
for(nm in names(columns)){tables[[nm]]<-read_table(nm,columns[[nm]],ids=ids);message(nm,': ',nrow(tables[[nm]]),' local rows')}
tables$code_status<-read_table('code_status',c('patient_id','start_dttm','code_status_category'),patients=patients)
rawg<-read_table('patient_assessments_raw_gcs',c('hospitalization_id','recorded_dttm','assessment_category','numerical_value'),ids=ids)
rawg$assessment_category[tolower(rawg$assessment_category)=='gcs_eyes']<-'gcs_eye'
rawg<-rawg[tolower(rawg$assessment_category)%in%c('gcs_eye','gcs_verbal','gcs_motor'),];rawg$assessment_category<-tolower(rawg$assessment_category)
rawg$numerical_value<-num(rawg$numerical_value)
valid<-with(rawg,(assessment_category=='gcs_eye' & numerical_value%in%1:4)|(assessment_category=='gcs_verbal'&numerical_value%in%1:5)|(assessment_category=='gcs_motor'&numerical_value%in%1:6))
g<-rawg[which(valid),] |> dplyr::group_by(hospitalization_id,recorded_dttm,assessment_category) |> dplyr::summarise(numerical_value=min(numerical_value),.groups='drop') |> tidyr::pivot_wider(names_from=assessment_category,values_from=numerical_value)
g<-with_columns(g,c('gcs_eye','gcs_verbal','gcs_motor'));g<-g[stats::complete.cases(g[c('gcs_eye','gcs_verbal','gcs_motor')]),]
gtot<-g |> dplyr::transmute(hospitalization_id,recorded_dttm,assessment_category='gcs_total',numerical_value=gcs_eye+gcs_verbal+gcs_motor)
imputed_gcs_rows<-sum(grepl('^gcs_',tables$patient_assessments$assessment_category))
tables$patient_assessments<-dplyr::bind_rows(tables$patient_assessments[!grepl('^gcs_',tables$patient_assessments$assessment_category),],rawg,gtot)
d<-d |> dplyr::left_join(sofa_features(d,tables,cfg),by='hospitalization_id')
land<-landmark_cohort(d,72)
flow<-land |> dplyr::count(landmark_reason,name='n');fed_write(flow,file.path(root,'aggregate/landmark_flow.csv'))
cohort_counts<-tibble::tibble(site='MIMIC',eligible_hospitalizations=nrow(d),patients=dplyr::n_distinct(d$patient_id),landmark_hospitalizations=sum(land$landmark_reason=='eligible'),previous_lung_cancer_hospitalizations=sum(d$prior_completed_lung_cancer_encounter),strict_poa_eligible=NA_integer_)
# Strict POA count is derived directly, not from the eligibility adapter.
cohort_counts$strict_poa_eligible<-length(unique(dx$hospitalization_id[which(is_lung & num(dx$poa_present)==1)]))
fed_write(cohort_counts,file.path(root,'aggregate/cohort_counts.csv'))
fed_json(list(site='MIMIC',source='MIMIC-IV-Ext-CLIF v1.1.0',cohort='Adult ICU hospitalizations with current-encounter C34.x/162.x primary-organ lung-cancer codes; POA unavailable',poa_nonmissing_rows=poa_known,poa_not_imputed=TRUE,eligibility_adapter_only=TRUE,calendar_window='all source periods; shifted dates not interpreted as real 2018-2024 calendar years',baseline_diagnosis_covariates='unavailable; discharge codes not substituted for baseline POA disease',gcs='Raw same-timestamp complete component triplets only; no best-score defaults or ETT-to-15 imputation',source_imputed_gcs_rows_excluded=imputed_gcs_rows,raw_complete_gcs_totals=nrow(gtot),primary_consensus='Frozen UCMC/NU/Rush full consensus; no MIMIC discovery or outcomes contribute'),file.path(root,'aggregate/compatibility_manifest.json'))
saveRDS(list(cohort=d,landmark=land,tables=tables,cfg=cfg),file.path(local,'prepared.rds'))
message('MIMIC coded-lung-cancer eligible=',nrow(d),'; +72h landmark=',sum(land$landmark_reason=='eligible'),'; prior completed lung-cancer encounter=',sum(d$prior_completed_lung_cancer_encounter))
