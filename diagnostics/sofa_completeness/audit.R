source('code/lib/load.R');source('diagnostics/sofa_completeness/recovery.R')
suppressPackageStartupMessages(library(arrow))
site<-commandArgs(trailingOnly=TRUE)[1]
source('diagnostics/sofa_completeness/context.R');ctx<-sofa_context(site);roots<-ctx$roots;paths<-ctx$paths;out<-ctx$output
a<-readRDS(file.path(roots[site],'local/prepared.rds'));d<-a$cohort;d<-d[d$hospitalization_id%in%readRDS(file.path(roots[site],'local/derivation.rds'))$analyses$h72_primary$data$hospitalization_id,]
t<-a$tables;cfg<-ctx$cfg;cfg$analysis<-a$parameters; ids<-d$hospitalization_id
window<-function(x,ts,pre=0){x |> dplyr::inner_join(d |> dplyr::select(hospitalization_id,t0,admission_dttm),by='hospitalization_id') |> dplyr::filter(.data[[ts]]>=pmax(admission_dttm,t0-pre*3600),.data[[ts]]<t0+24*3600)}
# Source categories/names/units are read locally to detect projected-away evidence.
labs<-sofa_source_rows(paths[site],'labs',ids)
labs<-with_columns(labs,c('lab_name','lab_loinc_code','reference_unit'))
for(ts in grep('_dttm$',names(labs),value=TRUE))labs[[ts]]<-utc(labs[[ts]])
labmeta<-labs |> dplyr::select(dplyr::any_of(c('lab_name','lab_category','reference_unit','lab_loinc_code'))) |> dplyr::distinct()
write.csv(labmeta,file.path(out,'reports',paste0(site,'_all_lab_metadata.csv')),row.names=FALSE)
nm<-toupper(labs$lab_name);lc<-as.character(labs$lab_loinc_code)
altpo<-labs$lab_category!='po2_arterial' & (lc%in%c('2703-7','19255-9','19254-2') | (grepl('PO2',nm) & grepl('ARTERIAL',nm))) & !grepl('CORD|VENOUS|SATUR',nm)
altbil<-labs$lab_category!='bilirubin_total' & lc=='1975-2' & !grepl('URINE|DIRECT|CONJUGATED|PERITONEAL',nm)
altpo[is.na(altpo)]<-FALSE;altbil[is.na(altbil)]<-FALSE
fixed<-labs;fixed$lab_category[altpo]<-'po2_arterial';fixed$lab_category[altbil]<-'bilirubin_total'
fiolab<-grepl('FIO2|FRACTION.*INSPIRED.*OXYGEN',nm) & !grepl('VENOUS|CORD',nm);fiolab[is.na(fiolab)]<-FALSE
f<-fixed[fiolab,,drop=FALSE];f$fio<-standard_fio2(f$lab_value_numeric)
f<-f[is.finite(f$fio),,drop=FALSE]
# FiO2 accompanies an arterial sample at exactly the same result timestamp.
po<-fixed[fixed$lab_category=='po2_arterial',c('hospitalization_id','lab_result_dttm'),drop=FALSE]
f<-f |> dplyr::inner_join(dplyr::distinct(po),by=c('hospitalization_id','lab_result_dttm'))
radd<-tibble::tibble(hospitalization_id=f$hospitalization_id,recorded_dttm=f$lab_result_dttm,device_category=NA_character_,fio2_set=f$fio,lpm_set=NA_real_)
# The synthetic FiO2 row must not replace device information. Pair arterial FiO2
# to the raw device separately, refusing ratios <200 without known support.
rec<-sofa_recover_tuples(t)
rv<-window(rec$map,'recorded_dttm');rg<-window(rec$gcs,'recorded_dttm')
no_map<-d$hospitalization_id[is.na(d$sofa_cv)];no_gcs<-d$hospitalization_id[is.na(d$sofa_cns)]
newt<-t
newt$vitals<-dplyr::bind_rows(t$vitals,rec$map[rec$map$hospitalization_id%in%no_map,])
newt$patient_assessments<-dplyr::bind_rows(t$patient_assessments,rec$gcs[rec$gcs$hospitalization_id%in%no_gcs,])
newt$labs<-fixed |> dplyr::select(dplyr::all_of(names(t$labs)))
# Preserve all original respiratory rows; lab FiO2 is evaluated separately below.
scores<-sofa_features(d,newt,cfg)
mx<-function(x)if(all(is.na(x)))NA_real_ else max(x,na.rm=TRUE)
reason<-list();extra<-list()
for(i in seq_len(nrow(d))){
 id<-d$hospitalization_id[i];start<-d$t0[i];end<-start+24*3600
 cut<-function(x,ts)x[which(x$hospitalization_id==id & x[[ts]]>=start & x[[ts]]<end),,drop=FALSE]
 l<-cut(fixed,'lab_result_dttm');r<-cut(t$respiratory_support,'recorded_dttm');v<-cut(t$vitals,'recorded_dttm');med<-cut(t$medication_admin_continuous,'admin_dttm');ff<-cut(radd,'recorded_dttm')
 pa<-l[l$lab_category=='po2_arterial' & is.finite(num(l$lab_value_numeric)) & num(l$lab_value_numeric)>=20 & num(l$lab_value_numeric)<=800,,drop=FALSE]
 pfscore<-numeric()
 if(nrow(pa)&&nrow(ff)){
  fi<-nearest_value(pa$lab_result_dttm,ff$recorded_dttm,ff$fio2_set,0)
  vent<-nearest_value(pa$lab_result_dttm,r$recorded_dttm,ifelse(support_state(r$device_category)=='UNKNOWN',NA_real_,as.numeric(support_state(r$device_category)%in%c('IMV','NIPPV','CPAP'))),1)==1
  pf<-num(pa$lab_value_numeric)/fi
  pfscore<-ifelse(pf<100 & vent,4,ifelse(pf<200 & vent,3,ifelse(pf<300,2,ifelse(pf<400,1,0))))
 }
 scores$sofa_resp[i]<-mx(c(scores$sofa_resp[i],pfscore))
 extra[[i]]<-data.frame(hospitalization_id=id,has_valid_pao2=nrow(pa)>0,has_raw_fio2=any(is.finite(standard_fio2(r$fio2_set))),has_lab_fio2=nrow(ff)>0,has_spo2=any(v$vital_category=='spo2' & is.finite(num(v$vital_value))),map_tuple=nrow(rv[rv$hospitalization_id==id,])>0,gcs_tuple=nrow(rg[rg$hospitalization_id==id,])>0)
}
dom<-c('sofa_resp','sofa_cv','sofa_coag','sofa_liver','sofa_renal','sofa_cns')
scores$sofa_observed_domains<-rowSums(!is.na(scores[,dom]));scores$sofa_total<-ifelse(scores$sofa_observed_domains==6,rowSums(scores[,dom]),NA_real_)
# Optional earlier-lab availability: window change is explicitly separate.
lw<-window(fixed,'lab_result_dttm',24)
summary<-dplyr::bind_rows(lapply(c(dom,'sofa_total'),function(c)data.frame(site=site,component=c,n=nrow(d),old_available=sum(!is.na(d[[c]])),recovered_available=sum(!is.na(scores[[c]])),newly_available=sum(is.na(d[[c]]) & !is.na(scores[[c]])),changed_existing=sum(!is.na(d[[c]]) & !is.na(scores[[c]]) & d[[c]]!=scores[[c]]))))
write.csv(summary,file.path(out,'reports',paste0(site,'_component_completeness.csv')),row.names=FALSE)
saveRDS(list(original=d,recovered=scores,evidence=dplyr::bind_rows(extra),parameters=list(window='[ICU,ICU+24h)',map='exact timestamp unambiguous SBP/DBP only when CV missing',gcs='exact timestamp complete E/V/M only when CNS missing',lab='explicit LOINC/name mappings',fio2='same-result-time arterial panel FiO2; no imputation'),source_metadata=labmeta),file.path(out,'local',paste0(site,'_audit.rds')))
write.csv(data.frame(site=site,evidence=names(extra[[1]])[-1],hospitalizations=colSums(dplyr::bind_rows(extra)[,-1])),file.path(out,'reports',paste0(site,'_evidence_availability.csv')),row.names=FALSE)
write.csv(lw |> dplyr::filter(lab_category%in%c('bilirubin_total','platelet_count','creatinine','po2_arterial')) |> dplyr::group_by(lab_category) |> dplyr::summarise(n=dplyr::n_distinct(hospitalization_id[is.finite(num(lab_value_numeric))]),.groups='drop') |> dplyr::mutate(site=site,window='24h pre-ICU through ICU+24h; admission bounded; availability only, not score'),file.path(out,'reports',paste0(site,'_expanded_lab_window.csv')),row.names=FALSE)
print(summary);cat('DONE',site,'\n')
