# Independently implemented APACHE-II-inspired APS+age; NOT standard APACHE II/IV.
APACHE_COMPONENTS<-c('temperature','map','heart_rate','respiratory_rate','oxygenation','acid_base','sodium','potassium','creatinine','hematology','wbc','gcs')
apache_adapted_spec<-function()list(version='clif-adapted-apache-II-v1',window_hours=24L,origin='first_ICU_arrival',score='12 APS components plus age; excludes unverified chronic-health points and acute-renal-failure doubling',hematology='all-site primary proxy: estimated hematocrit percent = 3 * measured hemoglobin g/dL; observed hematocrit separate availability score',missing='primary total requires all12 APS components and age; separate partial-sum plus12 missing-indicator sensitivity',oxygenation='same-specimen gases; prior FiO2 within1h or same-specimen lab FiO2; explicit room air0.21; high-FiO2 A-a at sea level/RQ0.8',acid_base='worst observed arterial pH; if absent, serum bicarbonate sensitivity substitution',gcs='valid observed total or simultaneous complete components; unsedated status unverified',full_APACHE_II=FALSE,APACHE_IV=FALSE)
apache_points<-function(x,component){
 x<-num(x);out<-rep(NA_real_,length(x));ok<-is.finite(x)
 score<-switch(component,
 temperature=function(v)ifelse(v>=41,4,ifelse(v>=39,3,ifelse(v>=38.5,1,ifelse(v>=36,0,ifelse(v>=34,1,ifelse(v>=32,2,ifelse(v>=30,3,4))))))),
 map=function(v)ifelse(v>=160,4,ifelse(v>=130,3,ifelse(v>=110,2,ifelse(v>=70,0,ifelse(v>=50,2,4))))),
 heart_rate=function(v)ifelse(v>=180,4,ifelse(v>=140,3,ifelse(v>=110,2,ifelse(v>=70,0,ifelse(v>=55,1,ifelse(v>=40,2,4)))))),
 respiratory_rate=function(v)ifelse(v>=50,4,ifelse(v>=35,3,ifelse(v>=25,1,ifelse(v>=12,0,ifelse(v>=10,1,ifelse(v>=6,2,4)))))),
 ph=function(v)ifelse(v>=7.7,4,ifelse(v>=7.6,3,ifelse(v>=7.5,1,ifelse(v>=7.33,0,ifelse(v>=7.25,2,ifelse(v>=7.15,3,4)))))),
 bicarbonate=function(v)ifelse(v>=52,4,ifelse(v>=41,3,ifelse(v>=32,1,ifelse(v>=22,0,ifelse(v>=18,2,ifelse(v>=15,3,4)))))),
 sodium=function(v)ifelse(v>=180,4,ifelse(v>=160,3,ifelse(v>=155,2,ifelse(v>=150,1,ifelse(v>=130,0,ifelse(v>=120,2,ifelse(v>=111,3,4))))))),
 potassium=function(v)ifelse(v>=7,4,ifelse(v>=6,3,ifelse(v>=5.5,1,ifelse(v>=3.5,0,ifelse(v>=3,1,ifelse(v>=2.5,2,4)))))),
 creatinine=function(v)ifelse(v>=3.5,4,ifelse(v>=2,3,ifelse(v>=1.5,2,ifelse(v>=.6,0,2)))),
 hematocrit=function(v)ifelse(v>=60,4,ifelse(v>=50,2,ifelse(v>=46,1,ifelse(v>=30,0,ifelse(v>=20,2,4))))),
 wbc=function(v)ifelse(v>=40,4,ifelse(v>=20,2,ifelse(v>=15,1,ifelse(v>=3,0,ifelse(v>=1,2,4))))),
 gcs=function(v)ifelse(v>=3 & v<=15 & v%%1==0,15-v,NA_real_),
 age=function(v)ifelse(v<18|v>120,NA_real_,ifelse(v>=75,6,ifelse(v>=65,5,ifelse(v>=55,3,ifelse(v>=45,2,0))))),
 stop('Unknown APACHE component'))
 out[ok]<-score(x[ok]);out
}
apache_oxygen_points<-function(po2,co2,fio2){
 n<-length(po2);stopifnot(length(co2)==n,length(fio2)==n)
 p<-rep(NA_real_,n);valid<-is.finite(po2)&po2>=20&po2<=800&is.finite(fio2)&fio2>=.21&fio2<=1
 low<-valid&fio2<.5;p[low]<-ifelse(po2[low]<55,4,ifelse(po2[low]<=60,3,ifelse(po2[low]<=70,1,0)))
 high<-valid&fio2>=.5&is.finite(co2)&co2>=10&co2<=200
 aa<-fio2[high]*713-co2[high]/.8-po2[high]
 p[high]<-ifelse(aa>=500,4,ifelse(aa>=350,3,ifelse(aa>=200,2,0)))
 p
}
apache_prior<-function(target,time,value,tolerance=3600){
 out<-rep(NA_real_,length(target));ok<-!is.na(time)&is.finite(value);time<-as.numeric(time[ok]);value<-value[ok]
 if(!length(time))return(out)
 z<-data.frame(time,value) |> dplyr::group_by(time) |> dplyr::summarise(n=dplyr::n_distinct(value),value=dplyr::first(value),.groups='drop') |> dplyr::filter(n==1) |> dplyr::arrange(time)
 if(!nrow(z))return(out)
 j<-findInterval(as.numeric(target),z$time);ix<-which(j>0 & !is.na(target));good<-ix[as.numeric(target[ix])-z$time[j[ix]]<=tolerance]
 out[good]<-z$value[j[good]];out
}
apache_lab_normalize<-function(l){
 l<-with_columns(l,c('lab_name','reference_unit','lab_loinc_code','lab_collect_dttm'))
 cat<-tolower(l$lab_category);cat[is.na(cat)]<-'unmapped';nm<-toupper(l$lab_name);loinc<-as.character(l$lab_loinc_code)
 arterial<-!grepl('CORD|VENOUS|SATUR',nm)
 ix<-((loinc%in%c('2703-7','19254-2','19255-9'))|(grepl('PO2',nm)&grepl('ARTERIAL',nm)))&arterial;ix[is.na(ix)]<-FALSE;cat[ix]<-'po2_arterial'
 ix<-grepl('FIO2|FRACTION.*INSPIRED.*OXYGEN',nm)&!grepl('CORD|VENOUS',nm);ix[is.na(ix)]<-FALSE;cat[ix]<-'fio2_lab'
 ix<-(loinc%in%c('4544-3','20570-8')|grepl('HEMATOCRIT|HAEMATOCRIT',nm))&!grepl('URINE|PLEURAL|PERITONEAL',nm);ix[is.na(ix)]<-FALSE;cat[ix]<-'hematocrit'
 val<-num(l$lab_value_numeric);u<-tolower(gsub('[[:space:]]','',gsub('[µμ]','u',l$reference_unit)))
 # UCUM mm[Hg] denotes the same pressure as mmHg; preserve rejection of other units.
 u[!is.na(u)&u=='mm[hg]']<-'mmhg'
 missing<-is.na(u)|!nzchar(u)
 accepted<-rep(TRUE,nrow(l));assumed<-missing
 for(k in c('hemoglobin','hematocrit','creatinine','po2_arterial','pco2_arterial','sodium','potassium','bicarbonate','wbc')){
  ix<-which(cat==k);if(!length(ix))next
  allowed<-switch(k,hemoglobin=c('g/dl','g/l'),hematocrit=c('%','percent','l/l'),creatinine=c('mg/dl','umol/l'),po2_arterial=c('mmhg','kpa'),pco2_arterial=c('mmhg','kpa'),sodium=c('mmol/l','meq/l'),potassium=c('mmol/l','meq/l'),bicarbonate=c('mmol/l','meq/l'),wbc=c('10^3/ul','10*3/ul','10^9/l','k/ul','x10^3/ul','x10(3)/ul','thou/ul'))
  accepted[ix]<-missing[ix]|u[ix]%in%allowed
 }
 val[cat=='hemoglobin'&u=='g/l'&!is.na(u)]<-val[cat=='hemoglobin'&u=='g/l'&!is.na(u)]/10
 val[cat=='hematocrit'&u=='l/l'&!is.na(u)]<-val[cat=='hematocrit'&u=='l/l'&!is.na(u)]*100
 val[cat=='creatinine'&u=='umol/l'&!is.na(u)]<-val[cat=='creatinine'&u=='umol/l'&!is.na(u)]/88.4
 ix<-cat%in%c('po2_arterial','pco2_arterial')&!is.na(u)&u=='kpa';val[ix]<-val[ix]*7.50062
 val[!accepted]<-NA_real_;val[cat=='fio2_lab']<-standard_fio2(val[cat=='fio2_lab'])
 l$apache_category<-cat;l$apache_value<-val;l$apache_unit_assumed<-assumed;l$apache_unit_rejected<-!accepted
 l
}
apache_adapted_features<-function(d,tables){
 source('diagnostics/sofa_completeness/recovery.R',local=TRUE)
 rec<-sofa_recover_tuples(tables)
 lbs<-apache_lab_normalize(tables$labs)
 by<-function(x)split(x,x$hospitalization_id)
 vs<-by(tables$vitals);ls<-by(lbs);rs<-by(tables$respiratory_support);gs<-by(tables$patient_assessments);maps<-by(rec$map);gcss<-by(rec$gcs)
 mx<-function(x)if(!length(x)||all(is.na(x)))NA_real_ else max(x,na.rm=TRUE)
 rows<-list()
 for(i in seq_len(nrow(d))){
  id<-d$hospitalization_id[i];start<-d$t0[i];end<-start+86400
  cut<-function(x,ts,empty){if(is.null(x))return(empty[0,]);x[which(!is.na(x[[ts]])&x[[ts]]>=start&x[[ts]]<end),,drop=FALSE]}
  v<-cut(vs[[id]],'recorded_dttm',tables$vitals);l<-cut(ls[[id]],'lab_result_dttm',lbs);r<-cut(rs[[id]],'recorded_dttm',tables$respiratory_support);g<-cut(gs[[id]],'recorded_dttm',tables$patient_assessments)
  vv<-function(k,lo,hi){x<-num(v$vital_value[tolower(v$vital_category)==k]);x[!is.finite(x)|x<lo|x>hi]<-NA_real_;x}
  lv<-function(k,lo,hi){x<-l$apache_value[l$apache_category==k];x[!is.finite(x)|x<lo|x>hi]<-NA_real_;x}
  map<-vv('map',20,250);map_derived<-FALSE
  if(!any(is.finite(map))){q<-cut(maps[[id]],'recorded_dttm',rec$map);map<-q$vital_value;map_derived<-any(is.finite(map))}
  gc<-num(g$numerical_value[tolower(g$assessment_category)=='gcs_total']);gc[!is.finite(gc)|gc<3|gc>15|gc%%1!=0]<-NA_real_;gcs_derived<-FALSE
  if(!any(is.finite(gc))){q<-cut(gcss[[id]],'recorded_dttm',rec$gcs);gc<-q$numerical_value;gcs_derived<-any(is.finite(gc))}
  ph<-lv('ph_arterial',6.5,8);use_bicarb<-!any(is.finite(ph));acid<-if(use_bicarb)apache_points(lv('bicarbonate',2,80),'bicarbonate') else apache_points(ph,'ph')
  # Sample identity uses collection timestamp when valid; otherwise result timestamp.
  l$sample_key<-as.numeric(l$lab_result_dttm)
  collect<-utc(l$lab_collect_dttm);ok<-!is.na(collect)&collect>=start&collect<end&collect<=l$lab_result_dttm
  l$sample_key[ok]<-as.numeric(collect[ok])
  pa<-l[l$apache_category=='po2_arterial'&is.finite(l$apache_value)&l$apache_value>=20&l$apache_value<=800,,drop=FALSE]
  oxygen<-numeric();paired<-0L
  if(nrow(pa)){
   fio<-standard_fio2(r$fio2_set);state<-support_state(r$device_category);fio[is.na(fio)&state=='ROOM_AIR']<-.21
   fi<-apache_prior(utc(pa$sample_key),r$recorded_dttm,fio)
   co<-rep(NA_real_,nrow(pa))
   for(j in seq_len(nrow(pa))){
    q<-l[l$sample_key==pa$sample_key[j],,drop=FALSE]
    f<-unique(q$apache_value[q$apache_category=='fio2_lab'&is.finite(q$apache_value)])
    if(length(f)==1)fi[j]<-f
    c<-unique(q$apache_value[q$apache_category=='pco2_arterial'&is.finite(q$apache_value)])
    if(length(c)==1)co[j]<-c
   }
   oxygen<-apache_oxygen_points(pa$apache_value,co,fi);paired<-sum(is.finite(oxygen))
  }
  hb<-lv('hemoglobin',1,25);hct<-lv('hematocrit',5,80)
  points<-c(temperature=mx(apache_points(vv('temp_c',25,45),'temperature')),map=mx(apache_points(map,'map')),heart_rate=mx(apache_points(vv('heart_rate',10,300),'heart_rate')),respiratory_rate=mx(apache_points(vv('respiratory_rate',1,100),'respiratory_rate')),oxygenation=mx(oxygen),acid_base=mx(acid),sodium=mx(apache_points(lv('sodium',80,220),'sodium')),potassium=mx(apache_points(lv('potassium',1,12),'potassium')),creatinine=mx(apache_points(lv('creatinine',.1,25),'creatinine')),hematology=mx(apache_points(3*hb,'hematocrit')),wbc=mx(apache_points(lv('wbc',0,500),'wbc')),gcs=mx(apache_points(gc,'gcs')))
  age<-apache_points(d$age_years[i],'age');observed<-sum(is.finite(points));complete<-observed==12&&is.finite(age)
  hpoints<-points;hpoints['hematology']<-mx(apache_points(hct,'hematocrit'))
  row<-tibble::tibble(hospitalization_id=id,apache_age_points=age,apache_observed_components=observed,apache_adapted_complete=complete,apache_adapted_total=if(complete)sum(points)+age else NA_real_,apache_partial_sum=sum(points,na.rm=TRUE)+age,apache_observed_hct_total=if(all(is.finite(hpoints))&&is.finite(age))sum(hpoints)+age else NA_real_,apache_hb_proxy_used=is.finite(points['hematology']),apache_map_derived=map_derived,apache_gcs_derived=gcs_derived,apache_bicarbonate_substituted=use_bicarb&&is.finite(points['acid_base']),apache_oxygen_pairs=paired,apache_rejected_unit_rows=sum(l$apache_unit_rejected),apache_assumed_unit_rows=sum(l$apache_unit_assumed),apache_chronic_points=NA_real_,apache_acute_renal_modifier=NA_real_)
  for(k in names(points))row[[paste0('apache_',k)]]<-points[k]
  rows[[i]]<-row
 }
 dplyr::bind_rows(rows)
}
# Enrich frozen numeric rows with source metadata only after an exact row-key match.
# Source release updates cannot introduce different measurements into the comparator.
apache_enrich_cached_labs<-function(cached,source){
 source<-with_columns(source,c('lab_name','lab_loinc_code','reference_unit','lab_collect_dttm'))
 for(ts in grep('_dttm$',names(source),value=TRUE))source[[ts]]<-utc(source[[ts]])
 key<-function(x)paste(x$hospitalization_id,format(as.numeric(x$lab_result_dttm),digits=17,scientific=TRUE,trim=TRUE),x$lab_category,format(num(x$lab_value_numeric),digits=17,scientific=TRUE,trim=TRUE),sep='|')
 rawkey<-key(source);cachekey<-key(cached)
 # Reject conflicting source metadata at the same frozen numeric key.
 source$.key<-rawkey;source$.row<-seq_len(nrow(source))
 normalized<-apache_lab_normalize(source);source$.type<-normalized$apache_category
 meta<-source |> dplyr::group_by(.key) |> dplyr::summarise(row=dplyr::first(.row),types=dplyr::n_distinct(.type),units=dplyr::n_distinct(reference_unit),.groups='drop')
 ix<-match(cachekey,meta$.key);ri<-meta$row[ix]
 out<-cached
 for(c in c('lab_name','lab_loinc_code','reference_unit','lab_collect_dttm'))out[[c]]<-source[[c]][ri]
 conflict<-!is.na(ix)&(meta$types[ix]>1|meta$units[ix]>1)
 out$lab_name[conflict]<-NA_character_;out$lab_loinc_code[conflict]<-NA_character_;out$reference_unit[conflict]<-'ambiguous_source_unit'
 list(labs=out,audit=data.frame(cached_rows=nrow(cached),metadata_matched_rows=sum(!is.na(ix)),metadata_unmatched_rows=sum(is.na(ix)),conflicting_metadata_rows=sum(conflict),source_new_numeric_rows_ignored=sum(!rawkey%in%cachekey)))
}
# Portable component capture used by every site's supplement phase.
capture_apache_adapted_components<-function(prepared,analyses,cfg,root){
 source('diagnostics/sofa_completeness/context.R',local=TRUE)
 d<-analyses$h72_primary$data;tables<-prepared$tables
 raw<-sofa_source_rows(cfg$tables_path,'labs',d$hospitalization_id)
 cached<-tables$labs[tables$labs$hospitalization_id%in%d$hospitalization_id,]
 enriched<-apache_enrich_cached_labs(cached,raw);tables$labs<-enriched$labs
 # Audit actual first-24h inputs, preserving source unit spellings and rejection flags.
 audit<-apache_lab_normalize(enriched$labs) |>
  dplyr::inner_join(d |> dplyr::select(hospitalization_id,t0),by='hospitalization_id') |>
  dplyr::filter(lab_result_dttm>=t0,lab_result_dttm<t0+86400) |>
  dplyr::count(lab_category,apache_category,reference_unit,apache_unit_rejected,apache_unit_assumed,name='n')
 supplement_write(audit,root,'apache_adapted_lab_units',cfg,'n',c('apache_category','reference_unit'))
 a<-apache_adapted_features(d,tables)
 dir.create(file.path(root,'local/apache_adapted'),recursive=TRUE,showWarnings=FALSE)
 saveRDS(a,file.path(root,'local/apache_adapted/components.rds'))
 readr::write_csv(enriched$audit,file.path(root,'local/apache_adapted/source_metadata_match_audit.csv'))
 rows<-dplyr::bind_rows(lapply(c(APACHE_COMPONENTS,'adapted_total','observed_hct_total'),function(k){x<-a[[paste0('apache_',k)]];tibble::tibble(component=k,n=nrow(a),observed_n=sum(is.finite(x)),missing_n=sum(!is.finite(x)),method='CLIF_adapted_APACHE_II_APS_age_not_validated_APACHE')}))
 supplement_write(rows,root,'apache_adapted_components',cfg,c('n','observed_n','missing_n'),'component')
 flags<-c('hb_proxy_used','map_derived','gcs_derived','bicarbonate_substituted')
 z<-dplyr::bind_rows(lapply(flags,function(k)tibble::tibble(flag=k,n=nrow(a),flag_n=sum(a[[paste0('apache_',k)]]),other_n=sum(!a[[paste0('apache_',k)]]))))
 supplement_write(z,root,'apache_adapted_provenance',cfg,c('n','flag_n','other_n'),'flag')
 limits<-tibble::tibble(item=c('score','chronic_health','acute_renal_doubling','hematology','missing_components','model_status'),status=c('modified_APS_plus_age_not_APACHE_II_or_IV','unverified_excluded_not_assumed_absent','unverified_not_applied','Hb_times3_proxy_consistent_across_sites_observed_Hct_separate','unknown_primary_requires_all12_plus_age','components_only_heldout_model_phase_separate'))
 supplement_write(limits,root,'apache_adapted_status',cfg)
 invisible(a)
}
