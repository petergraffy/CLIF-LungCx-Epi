source('code/lib/load.R');suppressPackageStartupMessages(library(arrow))
site<-commandArgs(trailingOnly=TRUE)[1]
source('diagnostics/sofa_completeness/context.R');ctx<-sofa_context(site);roots<-ctx$roots;paths<-ctx$paths;out<-ctx$output
a<-readRDS(file.path(roots[site],'local/prepared.rds'));d<-a$cohort
landmark_ids<-readRDS(file.path(roots[site],'local/derivation.rds'))$analyses$h72_primary$data$hospitalization_id;d<-d[d$hospitalization_id%in%landmark_ids,];ids<-d$hospitalization_id
x<-sofa_source_rows(paths[site],'patient_assessments',ids) |> filter(assessment_category%in%c('gcs_total','gcs_eye','gcs_verbal','gcs_motor'))
x$recorded_dttm<-utc(x$recorded_dttm);x<-x |> inner_join(d |> select(hospitalization_id,t0),by='hospitalization_id') |> filter(recorded_dttm>=t0,recorded_dttm<t0+24*3600)
x<-with_columns(x,c('categorical_value','text_value'));n<-num(x$numerical_value)
# Only exact integer strings are candidates, never extract digits from clinical text.
parse_exact<-function(v){v<-trimws(as.character(v));ifelse(grepl('^[0-9]{1,2}$',v),num(v),NA_real_)}
x$candidate<-dplyr::coalesce(n,parse_exact(x$categorical_value),parse_exact(x$text_value))
x$numeric_missing<-is.na(n)
summary<-x |> group_by(assessment_category) |> summarise(stays=n_distinct(hospitalization_id),missing_numeric_rows=sum(numeric_missing),recoverable_string_rows=sum(numeric_missing & is.finite(candidate)),string_recoverable_stays=n_distinct(hospitalization_id[numeric_missing & is.finite(candidate)]),.groups='drop') |> mutate(site=site)
# Weight/medication and oxygenation audits within exactly the same ICU window.
win<-function(z,ts)z |> inner_join(d |> select(hospitalization_id,t0,admission_dttm),by='hospitalization_id') |> filter(.data[[ts]]>=t0,.data[[ts]]<t0+24*3600)
v<-win(a$tables$vitals,'recorded_dttm');m<-win(a$tables$medication_admin_continuous,'admin_dttm');l<-win(a$tables$labs,'lab_result_dttm');r<-win(a$tables$respiratory_support,'recorded_dttm')
units<-m |> filter(med_category%in%c('norepinephrine','epinephrine','dopamine','dobutamine'),num(med_dose)>0) |> group_by(med_category,med_dose_unit) |> summarise(rows=n(),stays=n_distinct(hospitalization_id),.groups='drop') |> mutate(site=site)
rawmed<-sofa_source_rows(paths[site],'medication_admin_continuous',ids);rawmed<-with_columns(rawmed,'med_name');rawmed$admin_dttm<-utc(rawmed$admin_dttm);rawmed<-win(rawmed,'admin_dttm')
nm<-toupper(rawmed$med_name)
rawmed$expected<-case_when(grepl('NOREPINEPHRINE',nm)~'norepinephrine',grepl('EPINEPHRINE',nm)~'epinephrine',grepl('DOBUTAMINE',nm)~'dobutamine',grepl('DOPAMINE',nm)~'dopamine',TRUE~NA_character_)
mismatch<-rawmed |> filter(!is.na(expected),expected!=tolower(med_category),num(med_dose)>0) |> group_by(med_name,med_category,expected,med_dose_unit) |> summarise(rows=n(),stays=n_distinct(hospitalization_id),.groups='drop') |> mutate(site=site)
# Hospitalization-level evidence stays local, summary counts have no identifiers.
weights<-a$tables$vitals |> filter(vital_category=='weight_kg',num(vital_value)>10,num(vital_value)<500) |> inner_join(d |> select(hospitalization_id,t0,admission_dttm),by='hospitalization_id') |> filter(recorded_dttm>=pmax(admission_dttm,t0-24*3600),recorded_dttm<t0+24*3600)
counts<-data.frame(site=site,unconverted_cv=sum(d$unconverted_cv_dose),pre24_or_first24_weight_stays=n_distinct(weights$hospitalization_id),first24_weight_stays=n_distinct(v$hospitalization_id[v$vital_category=='weight_kg' & num(v$vital_value)>10 & num(v$vital_value)<500]),med_mapping_mismatch_stays=n_distinct(rawmed$hospitalization_id[!is.na(rawmed$expected) & rawmed$expected!=tolower(rawmed$med_category) & num(rawmed$med_dose)>0]))
write.csv(summary,file.path(out,'reports',paste0(site,'_gcs_value_audit.csv')),row.names=FALSE);write.csv(units,file.path(out,'reports',paste0(site,'_vasoactive_units.csv')),row.names=FALSE);write.csv(mismatch,file.path(out,'reports',paste0(site,'_vasoactive_mapping_flags.csv')),row.names=FALSE);write.csv(counts,file.path(out,'reports',paste0(site,'_extra_source_counts.csv')),row.names=FALSE)
saveRDS(list(gcs=x,weights=weights,rawmed=rawmed),file.path(out,'local',paste0(site,'_extra_source.rds')));print(summary);print(counts)
