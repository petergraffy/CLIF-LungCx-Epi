source('diagnostics/revision20261022/lib.R')
site<-commandArgs(trailingOnly=TRUE)[1];home<-file.path(rev_root(),'sites',site);out<-file.path(home,'export');dir.create(file.path(home,'local'),recursive=TRUE,showWarnings=FALSE)
a<-readRDS(file.path(old_root(),'sites',site,'local/discovery.rds'))
d<-readRDS(file.path(old_root(),'sites',site,'local/full/outcomes.rds'))$data
p<-a$panel |> dplyr::arrange(hospitalization_id,h)
# Hourly evidence states remain distinct: H1_CU is not hypoxemia-only.
p$arf_type<-dplyr::case_when(p$hypox=='H1' & p$hyper=='C1'~'mixed',p$hypox=='H1' & p$hyper=='C0'~'hypoxemia_with_negative_hypercapnia',p$hypox=='H0' & p$hyper=='C1'~'hypercapnia_with_negative_hypoxemia',p$hypox=='H1'~'hypoxemia_other_channel_unknown',p$hyper=='C1'~'hypercapnia_other_channel_unknown',p$hypox=='H0' & p$hyper=='C0'~'both_negative',TRUE~'uncertain_without_positive_evidence')
arf<-p |> dplyr::group_by(hospitalization_id) |> dplyr::summarise(hypox_fraction=mean(hypox=='H1'),hyper_fraction=mean(hyper=='C1'),mixed_fraction=mean(hypox=='H1' & hyper=='C1'),any_hypox=any(hypox=='H1'),any_hyper=any(hyper=='C1'),any_mixed=any(hypox=='H1' & hyper=='C1'),first_detected_hour=if(any(hypox=='H1'|hyper=='C1'))min(h[hypox=='H1'|hyper=='C1']) else NA_real_,positive_fraction=mean(hypox=='H1'|hyper=='C1'),late_positive=any((hypox=='H1'|hyper=='C1') & h>=48),.groups='drop')
d<-d |> dplyr::left_join(arf,by='hospitalization_id')
d$arf_pattern<-with(d,ifelse(any_mixed,'simultaneous_mixed',ifelse(any_hypox & any_hyper,'both_detected_at_different_hours',ifelse(any_hypox,'hypoxemia_detected_only',ifelse(any_hyper,'hypercapnia_detected_only','no_positive_evidence')))))
d$arf_timing<-ifelse(is.na(d$first_detected_hour),'no_positive_evidence',ifelse(d$first_detected_hour<24,'first_0_24h',ifelse(d$first_detected_hour<48,'first_24_48h','first_48_72h')))
d$end_code<-ifelse(d$code_at_landmark=='unknown','unknown',ifelse(d$code_at_landmark=='full','full','nonfull'))
bands<-list()
for(v in c('code_at_admission','code_at_landmark','code_changed','new_dni','new_comfort','arf_pattern','arf_timing','late_positive','metastatic_cancer','malignant_pleural_effusion','cachexia','covid','any_imv','any_vasopressor'))bands[[v]]<-band_table(d,'phenotype',v,d[[v]])
trans<-band_table(d,c('phenotype','code_at_admission'),'admission_to_landmark_code',d$code_at_landmark)
fed_write(dplyr::bind_rows(bands) |> dplyr::mutate(site=site),file.path(out,'phenotype_characteristics.csv'))
fed_write(trans |> dplyr::mutate(site=site),file.path(out,'code_transitions.csv'))
# Each type/hour category contributes at most one count per hospitalization/bin.
z<-p |> dplyr::left_join(d |> dplyr::select(hospitalization_id,phenotype),by='hospitalization_id') |> dplyr::mutate(bin=h%/%6) |> dplyr::group_by(hospitalization_id,phenotype,bin) |> dplyr::summarise(category=fed_mode(arf_type),.groups='drop')
fed_write(band_table(z,c('phenotype','bin'),'arf_evidence_type',z$category) |> dplyr::mutate(site=site),file.path(out,'arf_type_profiles.csv'))
# Local clinical table has explicit denominators. Export only n/event bands.
continuous<-d |> dplyr::group_by(phenotype) |> dplyr::summarise(n=dplyr::n(),dplyr::across(c(age_years,charlson_noncancer,sofa_total,hours_IMV,hours_HFNC,hours_NIPPV,hours_UNKNOWN,hypox_fraction,hyper_fraction,mixed_fraction,positive_fraction),~mean(.x,na.rm=TRUE)),.groups='drop')
for(i in seq_len(nrow(continuous)))for(v in setdiff(names(continuous),c('phenotype','n')))if(sum(!is.na(d[[v]][d$phenotype==continuous$phenotype[i]]))<aggregate_min_cell(10))continuous[[v]][i]<-NA_real_
bands<-t(vapply(continuous$n,fed_band_n,numeric(2)));continuous$n_lower<-bands[,1];continuous$n_upper<-bands[,2]
for(i in which(continuous$n<aggregate_min_cell(10)))continuous[i,setdiff(names(continuous),c('phenotype','n','n_lower','n_upper'))]<-NA_real_
continuous$n<-NULL;fed_write(continuous |> dplyr::mutate(site=site),file.path(out,'clinical_table.csv'))
# Actual model denominators are retained locally, with banded public metadata.
subsets<-list(all=rep(TRUE,nrow(d)),full_at_admission=d$code_at_admission=='full',full_through_landmark=d$code_at_admission=='full' & d$code_at_landmark=='full' & !is.na(d$code_changed) & !d$code_changed & !is.na(d$new_dni) & !d$new_dni & !is.na(d$new_comfort) & !d$new_comfort,no_covid=!d$covid,first_hospitalization=d$first_hospitalization,arf_or_high_support=d$measured_arf|d$high_support)
rows<-list();meta<-list();fits<-list();exact<-list()
for(sub in names(subsets))for(spec in if(sub=='all')c('core','end_code','arf_channels') else 'core')for(endpoint in cox_endpoints()){
 z<-d[which(subsets[[sub]]),];z$event<-cox_event(z,endpoint)
 covars<-model_spec()$core
 if(spec=='end_code'){
 z$landmark_nonfull<-as.numeric(z$end_code=='nonfull');z$landmark_unknown<-as.numeric(z$end_code=='unknown');z$new_dni_observed<-as.numeric(!is.na(z$new_dni)&z$new_dni);z$new_comfort_observed<-as.numeric(!is.na(z$new_comfort)&z$new_comfort);covars<-c(covars,'landmark_nonfull','landmark_unknown','new_dni_observed','new_comfort_observed')
 }
 if(spec=='arf_channels')covars<-c(covars,'hypox_fraction','hyper_fraction','mixed_fraction')
 terms<-c(paste0('phenotype_',setdiff(revision_phenotypes(),'P1')),covars);terms<-terms[vapply(z[terms],function(x)length(unique(x))>1,logical(1))]
 f<-as.formula(paste('survival::Surv(followup,event)~',paste(c(terms,'strata(info_status)','cluster(patient_id)'),collapse='+')))
 warn<-character();fit<-tryCatch(withCallingHandlers(survival::coxph(f,data=z),warning=function(w){warn<<-c(warn,conditionMessage(w));invokeRestart('muffleWarning')}),error=identity)
 good<-inherits(fit,'coxph') && !length(warn) && all(is.finite(coef(fit)))
 key<-paste(sub,spec,endpoint,sep='_');fits[[key]]<-fit
 exact[[key]]<-tibble::tibble(site=site,subset=sub,model=spec,endpoint=endpoint,n=nrow(z),patients=dplyr::n_distinct(z$patient_id),events=sum(z$event),status=if(good)'ok' else paste(c(if(inherits(fit,'error'))conditionMessage(fit),warn),collapse=';'))
 nb<-fed_band_n(nrow(z));eb<-fed_band_n(sum(z$event));meta[[key]]<-tibble::tibble(site=site,subset=sub,model=spec,endpoint=endpoint,n_lower=nb[1],n_upper=nb[2],events_lower=eb[1],events_upper=eb[2],status=if(good)'ok' else 'nonestimable',formula=paste(deparse(f),collapse=' '))
 for(ph in setdiff(revision_phenotypes(),'P1')){
 term<-paste0('phenotype_',ph);safe<-all(vapply(c('P1',ph),function(q){zz<-z[z$phenotype==q,];sum(zz$event)>=aggregate_min_cell(10) & sum(1-zz$event)>=aggregate_min_cell(10)},logical(1)))
 ok<-good && safe && term%in%names(coef(fit));rows[[paste(key,ph)]]<-tibble::tibble(site=site,subset=sub,model=spec,endpoint=endpoint,phenotype=ph,reference='P1',estimate=if(ok)unname(coef(fit)[term]) else NA_real_,std_error=if(ok)sqrt(vcov(fit)[term,term]) else NA_real_,status=if(ok)'ok' else if(!safe)'withheld_small_partition' else 'nonestimable')
 }
}
fed_write(dplyr::bind_rows(rows),file.path(out,'sensitivity_associations.csv'));fed_write(dplyr::bind_rows(meta),file.path(out,'sensitivity_metadata.csv'))
readr::write_csv(dplyr::bind_rows(exact),file.path(home,'local/model_denominators.csv'));saveRDS(list(data=d,models=fits),file.path(home,'local/review.rds'))
# Early events are described by the absorbing sensitivity's LOCAL clusters only.
wf<-readRDS(file.path(a$input_root,'local/waterfall_derivation.rds'))$analyses$h72_absorbing
ed<-wf$data;ed$early_event<-ifelse(ed$event_time<=ed$t0+72*3600,ed$event_type,'alive_in_hospital_at_72h')
fed_write(band_table(ed,'local_cluster','early_event',ed$early_event) |> dplyr::mutate(site=site,analysis='absorbing_descriptive_only'),file.path(out,'early_events.csv'))
# Reconcile source dispositions including invalid outcomes, never guess a death.
m<-jsonlite::fromJSON(file.path(a$input_root,'local/input_manifest.json'));prepared<-readRDS(file.path(m$baseline_root,'local/prepared.rds'))
co<-prepared$cohort;co$flow<-ifelse(!co$valid_outcome,'unclassifiable_outcome',ifelse(co$event_time<=co$t0+72*3600,paste0('early_',co$event_type),'primary_landmark'))
fed_write(band_table(co,character(),'cohort_flow',co$flow) |> dplyr::mutate(site=site),file.path(out,'cohort_flow.csv'))
# APACHE feasibility: measured components in first ICU 24h, without imputation.
t<-prepared$tables
available<-list();requirements<-list(temperature=c('vitals','temp_c'),map=c('vitals','map'),heart_rate=c('vitals','heart_rate'),respiratory_rate=c('vitals','respiratory_rate'),sodium=c('labs','sodium'),potassium=c('labs','potassium'),creatinine=c('labs','creatinine'),hematocrit=c('labs','hematocrit'),wbc=c('labs','wbc'),arterial_ph=c('labs','ph_arterial'),arterial_po2=c('labs','po2_arterial'),arterial_pco2=c('labs','pco2_arterial'),gcs=c('patient_assessments','gcs_total'))
for(nm in names(requirements)){
 r<-requirements[[nm]];x<-t[[r[1]]];catcol<-switch(r[1],vitals='vital_category',labs='lab_category',patient_assessments='assessment_category');ts<-if(r[1]=='labs')'lab_result_dttm' else 'recorded_dttm';val<-switch(r[1],vitals='vital_value',labs='lab_value_numeric',patient_assessments='numerical_value')
 x<-x[tolower(x[[catcol]])==r[2] & is.finite(num(x[[val]])),] |> dplyr::inner_join(d |> dplyr::select(hospitalization_id,t0),by='hospitalization_id')
 ids<-unique(x$hospitalization_id[x[[ts]]>=x$t0 & x[[ts]]<x$t0+24*3600]);n<-sum(d$hospitalization_id%in%ids);b<-fed_band_n(n)
 available[[nm]]<-tibble::tibble(site=site,component=nm,available_lower=b[1],available_upper=b[2],status=if(n==0)'absent' else 'available_in_some_hospitalizations')
}
for(nm in c('validated_apache_ii','validated_apache_iv','apache_chronic_health_criteria','apache_surgical_admission_class','apache_acute_renal_failure_modifier','apache_iv_admission_diagnosis_and_required_history'))available[[nm]]<-tibble::tibble(site=site,component=nm,available_lower=0,available_upper=0,status='not_available_in_current_extract')
fed_write(dplyr::bind_rows(available),file.path(out,'severity_availability.csv'))
fed_json(list(site=site,primary_n=nrow(d),primary_all_retained=TRUE,apache_ii_computable=FALSE,apache_iv_computable=FALSE,apache_reason='No validated score; chronic health, surgical context and acute renal failure modifier unverified. Hematocrit absent at UCMC. Do not label an incomplete APS proxy APACHE.',sofa='Existing first-24-hour measured-component modified SOFA; creatinine only renal domain; no missing-as-normal primary score',source_hashes=as.list(tools::md5sum(c('diagnostics/revision20261022/site_review.R','diagnostics/revision20261022/lib.R')))),file.path(out,'review_manifest.json'))
message(site,' clinical review and adjusted sensitivities complete')
