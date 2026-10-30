source('code/lib/load.R');source('diagnostics/sofa_completeness/recovery.R')
site<-commandArgs(trailingOnly=TRUE)[1]
source('diagnostics/sofa_completeness/context.R');ctx<-sofa_context(site);roots<-ctx$roots;out<-ctx$output
a<-readRDS(file.path(roots[site],'local/prepared.rds'));b<-readRDS(file.path(out,'local',paste0(site,'_audit.rds')));e<-readRDS(file.path(out,'local',paste0(site,'_extra_source.rds')))
d<-b$original;z<-b$recovered;dom<-c('sofa_resp','sofa_cv','sofa_coag','sofa_liver','sofa_renal','sofa_cns');ids<-d$hospitalization_id
# Work per hospitalization with raw, pre-waterfall observations only.
rby<-split(a$tables$respiratory_support,a$tables$respiratory_support$hospitalization_id);vby<-split(a$tables$vitals,a$tables$vitals$hospitalization_id);lby<-split(a$tables$labs,a$tables$labs$hospitalization_id);mby<-split(e$rawmed,e$rawmed$hospitalization_id);wby<-split(e$weights,e$weights$hospitalization_id)
mp<-sofa_recover_tuples(a$tables)$map;mpby<-split(mp,mp$hospitalization_id)
mx<-function(x)if(!length(x)||all(is.na(x)))NA_real_ else max(x,na.rm=TRUE)
strict_cv<-extended_cv<-sf_resp<-rep(NA_real_,nrow(d)); flags<-list()
for(i in seq_len(nrow(d))){
 id<-ids[i];start<-d$t0[i];end<-start+86400
 cut<-function(x,ts){if(is.null(x))return(NULL);x[which(x[[ts]]>=start & x[[ts]]<end),,drop=FALSE]}
 v<-cut(vby[[id]],'recorded_dttm');r<-cut(rby[[id]],'recorded_dttm');med<-mby[[id]];map<-cut(mpby[[id]],'recorded_dttm');lw<-lby[[id]]
 if(is.null(v))v<-a$tables$vitals[0,]
 if(is.null(med))med<-e$rawmed[0,]
 category<-tolower(med$med_category);validname<-!is.na(med$expected) & category%in%c('dextrose_5_water','dextrose')
 category[validname]<-med$expected[validname]
 dose_unit<-tolower(gsub(' ','',med$med_dose_unit));dose<-num(med$med_dose)
 # Use a measured first-24h weight first; sensitivity adds latest pre-ICU weight.
 w<-num(v$vital_value[v$vital_category=='weight_kg']);w<-w[is.finite(w)&w>10&w<500];weight<-if(length(w))median(w) else NA_real_
 prew<-wby[[id]];extended_weight<-weight
 if(is.na(extended_weight) && !is.null(prew)){
  prew<-prew[which(prew$recorded_dttm<start),,drop=FALSE]
  if(nrow(prew)){prew<-prew[order(prew$recorded_dttm,decreasing=TRUE),];extended_weight<-num(prew$vital_value[1])}
 }
 maps<-num(v$vital_value[v$vital_category=='map']);maps<-maps[is.finite(maps)&maps>=20&maps<=250]
 if(!length(maps) && !is.null(map))maps<-map$vital_value
 cvcalc<-function(we){
  converted<-dose_per_kg_min(dose,dose_unit,rep(we,length(dose)))
  # Dobutamine uses any dose, so weight conversion is not required for its score.
  dobutamine_known<-category=='dobutamine' & dose>0 & dose_unit%in%c('mcg/min','mcg/kg/min','mcg/kg/hr','mcg/kg/hour','mg/hr','mg/hour','mg/kg/hr','mg/kg/hour')
  unresolved<-any(dose>0 & category%in%c('dopamine','epinephrine','norepinephrine') & is.na(converted),na.rm=TRUE) || any(dose>0 & category=='dobutamine' & !dobutamine_known,na.rm=TRUE)
  dop<-mx(converted[category=='dopamine']);ep<-mx(converted[category=='epinephrine']);ne<-mx(converted[category=='norepinephrine'])
  if(unresolved)return(NA_real_)
  mx(c(ifelse(maps<70,1,0),ifelse(dop>15,4,ifelse(dop>5,3,ifelse(dop>0,2,0))),ifelse(ep>.1,4,ifelse(ep>0,3,0)),ifelse(ne>.1,4,ifelse(ne>0,3,0)),if(any(dobutamine_known,na.rm=TRUE))2 else NA_real_))
 }
 strict_cv[i]<-cvcalc(weight);extended_cv[i]<-cvcalc(extended_weight)
 sp<-v[v$vital_category=='spo2' & is.finite(num(v$vital_value)) & num(v$vital_value)>=50 & num(v$vital_value)<97,,drop=FALSE]
 if(nrow(sp) && !is.null(r) && nrow(r)){
  fio<-standard_fio2(r$fio2_set);state<-support_state(r$device_category);fio[is.na(fio)&state=='ROOM_AIR']<-.21
  fi<-nearest_value(sp$recorded_dttm,r$recorded_dttm,fio,1)
  vent<-nearest_value(sp$recorded_dttm,r$recorded_dttm,ifelse(state=='UNKNOWN',NA_real_,as.numeric(state%in%c('IMV','NIPPV','CPAP'))),1)==1
  spnum<-num(sp$vital_value);q<-11700/((100/spnum)-1);pao2<-(sqrt(q^2+50^3)+q)^(1/3)-(sqrt(q^2+50^3)-q)^(1/3);pf<-pao2/fi
  sf_resp[i]<-mx(ifelse(pf<100 & vent,4,ifelse(pf<200 & vent,3,ifelse(pf<300,2,ifelse(pf<400,1,0)))))
 }
 flags[[i]]<-data.frame(hospitalization_id=id,reclassified_drug_rows=sum(validname & dose>0,na.rm=TRUE),weight_recovered=is.na(weight)&is.finite(extended_weight),resp_no_abg=!b$evidence$has_valid_pao2[i],resp_no_score_despite_abg=b$evidence$has_valid_pao2[i]&is.na(z$sofa_resp[i]))
}
scoretotal<-function(x){x$sofa_observed_domains<-rowSums(!is.na(x[,dom]));x$sofa_total<-ifelse(x$sofa_observed_domains==6,rowSums(x[,dom]),NA_real_);x}
strict<-z;strict$sofa_cv<-strict_cv;strict<-scoretotal(strict)
weight_sens<-strict;weight_sens$sofa_cv<-extended_cv;weight_sens<-scoretotal(weight_sens)
spo2_sens<-weight_sens;use<-is.na(spo2_sens$sofa_resp);spo2_sens$sofa_resp[use]<-sf_resp[use];spo2_sens<-scoretotal(spo2_sens)
# This estimate is an explicit sensitivity, not observed PaO2 or standard SOFA.
variants<-list(original=d,measurement_recovery=z,reviewed_drug_mapping=strict,preICU_measured_weight_sensitivity=weight_sens,SpO2_estimated_PaO2_sensitivity=spo2_sens)
sumtab<-bind_rows(lapply(names(variants),function(v){x<-variants[[v]];bind_rows(lapply(c(dom,'sofa_total'),function(c)data.frame(site=site,variant=v,component=c,n=nrow(x),available=sum(!is.na(x[[c]])),changed_existing=sum(!is.na(d[[c]])&!is.na(x[[c]])&d[[c]]!=x[[c]]))))}))
write.csv(sumtab,file.path(out,'reports',paste0(site,'_staged_completeness.csv')),row.names=FALSE)
saveRDS(list(variants=variants,flags=bind_rows(flags),sf_estimate=sf_resp),file.path(out,'local',paste0(site,'_extended.rds')))
write.csv(data.frame(site=site,evidence=names(flags[[1]])[-1],count=colSums(bind_rows(flags)[,-1])),file.path(out,'reports',paste0(site,'_recovery_flags.csv')),row.names=FALSE)
print(sumtab[sumtab$component%in%c('sofa_cv','sofa_resp','sofa_total'),]);cat('DONE',site,'\n')
