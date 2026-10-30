# Publication descriptors; aggregate sufficient statistics, never patient rows.
descriptive_registry<-function(){
 cont<-c(age_years='Age (years)',charlson_noncancer='Noncancer Charlson index',sofa_total='Original modified SOFA',sofa_observed_domains='Observed SOFA domains',apache_adapted_total='CLIF-adapted APACHE APS+age proxy',ward_to_icu_hours='Hospital admission to first ICU (hours)',hospital_los_days='Hospital LOS (days; whole stay descriptive)',hours_ROOM_AIR='Room-air hours in first72h',hours_LOW_O2='Low-flow oxygen hours in first72h',hours_HFNC='HFNC hours in first72h',hours_CPAP='CPAP hours in first72h',hours_NIPPV='NIPPV hours in first72h',hours_IMV='IMV hours in first72h',hours_OTHER='Other support hours in first72h',hours_UNKNOWN='Unknown support hours in first72h',vasopressor_documented_hours='Documented vasopressor hours in first72h',support_coverage='Known-support fraction in first72h',physiology_unknown_hours='Unknown-physiology hours in first72h',hypox_positive_hours='Hypoxemia-positive hours in first72h',hyper_positive_hours='Hypercapnia-positive hours in first72h',mixed_positive_hours='Simultaneous-positive hours in first72h',hypox_unknown_hours='Hypoxemia-unknown hours in first72h',hyper_unknown_hours='Hypercapnia-unknown hours in first72h')
 cat<-c(sex_group='Sex (harmonized)',race_group='Race (White/Black/other/unknown)',ethnicity_group='Ethnicity (harmonized)',sex_category='Sex (source category)',race_category='Race (source category)',ethnicity_category='Ethnicity (source category)',admission_year='Hospital admission year',first_hospitalization='First eligible hospitalization at this site',advanced_cancer='POA advanced-cancer proxy',metastatic_cancer='POA metastatic-cancer code',malignant_pleural_effusion='POA malignant-pleural-effusion code',cachexia='POA cachexia code',covid='POA U07.1 code',code_at_admission='Code status at ICU admission',code_at_landmark='Code status at72h',code_changed='Code-status change through72h',new_dni='New documented DNI through72h',new_comfort='New documented comfort order through72h',first_hour_support='Support in first ICU hour',support_at_landmark='Support in last completed landmark hour',any_low_o2='Any low-flow oxygen through72h',any_hfnc='Any HFNC through72h',any_cpap='Any CPAP through72h',any_nippv='Any NIPPV through72h',any_imv='Any IMV through72h',any_vasopressor='Any documented vasopressor through72h',measured_arf='Any positive measured ARF evidence through72h',any_hypox='Any positive hypoxemia evidence',any_hyper='Any positive hypercapnia evidence',any_mixed='Any simultaneous positive evidence',whole_stay_outcome='Terminal hospital outcome (whole stay descriptive)',whole_stay_death='Hospital death (whole stay descriptive)',whole_stay_hospice='Hospice discharge (whole stay descriptive)',whole_stay_composite='Death or hospice (whole stay descriptive)',outcome28='Hospital event within28 days AFTER landmark',death28='Hospital death within28 days AFTER landmark',hospice28='Hospice within28 days AFTER landmark',composite28='Death or hospice within28 days AFTER landmark',histology='Histology (timing unverified)',stage='Cancer stage (timing unverified)',ecog='ECOG (timing unverified)',therapy_line='Therapy line (timing unverified)',therapy_intent='Therapy intent (timing unverified)',molecular_subtype='Molecular subtype (timing unverified)',treatment_type='Treatment type (timing unverified)')
 domains<-c(mi='MI',chf='Heart failure',pvd='Peripheral vascular disease',cevd='Cerebrovascular disease',dementia='Dementia',cpd='Chronic pulmonary disease',rheumd='Rheumatologic disease',pud='Peptic ulcer disease',mld='Mild liver disease',diab='Diabetes without complication',diabwc='Diabetes with complication',hp='Hemiplegia/paraplegia',rend='Renal disease',msld='Moderate/severe liver disease',aids='AIDS/HIV')
 cat<-c(cat,setNames(paste0('POA coded Charlson: ',domains),paste0('charlson_',names(domains))))
 z<-dplyr::bind_rows(tibble::tibble(variable=names(cont),label=unname(cont),type='continuous'),tibble::tibble(variable=names(cat),label=unname(cat),type='categorical'))
 z$primary_only<-z$variable%in%c('apache_adapted_total','code_at_admission','code_at_landmark','code_changed','new_dni','new_comfort','first_hour_support','support_at_landmark','any_low_o2','any_hfnc','any_cpap','any_nippv','any_imv','any_vasopressor','measured_arf','any_hypox','any_hyper','any_mixed','outcome28','death28','hospice28','composite28','vasopressor_documented_hours','support_coverage','physiology_unknown_hours','hypox_positive_hours','hyper_positive_hours','mixed_positive_hours','hypox_unknown_hours','hyper_unknown_hours')|grepl('^hours_',z$variable)
 z$note<-ifelse(grepl('^charlson_',z$variable)&z$variable!='charlson_noncancer','POA coding: false means no corresponding code, not confirmed absence of disease',ifelse(grepl('whole_stay|hospital_los',z$variable),'Descriptive only; includes death/discharge; not an inferential LOS endpoint',ifelse(z$variable=='apache_adapted_total','Modified proxy; requires all12 adapted components; not standard APACHE II/IV','')))
 z
}
descriptive_comorbidities<-function(co,dx){
 if(is.null(dx)||!all(c('poa_present','diagnosis_code','diagnosis_code_format')%in%names(dx)))return(tibble::tibble(hospitalization_id=co$hospitalization_id))
 dx<-dx |> dplyr::filter(num(poa_present)==1) |> dplyr::mutate(code=toupper(gsub('[^A-Za-z0-9]','',diagnosis_code)),version=dplyr::case_when(grepl('10',diagnosis_code_format)~'10',grepl('9',diagnosis_code_format)~'9',TRUE~NA_character_))
 rows<-list();for(v in c('10','9')){z<-dx |> dplyr::filter(version==v,!is.na(code),nzchar(code)) |> dplyr::select(id=hospitalization_id,code);if(nrow(z))rows[[v]]<-comorbidity::comorbidity(z,id='id',code='code',map=paste0('charlson_icd',v,'_quan'),assign0=TRUE)}
 if(!length(rows))return(tibble::tibble(hospitalization_id=co$hospitalization_id))
 z<-dplyr::bind_rows(rows);domains<-setdiff(names(z),'id')
 z<-z |> dplyr::group_by(id) |> dplyr::summarise(dplyr::across(dplyr::all_of(domains),max),.groups='drop')
 z$mld[z$msld==1]<-0;z$diab[z$diabwc==1]<-0;z$canc<-0;z$metacanc<-0
 weights<-c(mi=1,chf=1,pvd=1,cevd=1,dementia=1,cpd=1,rheumd=1,pud=1,mld=1,diab=1,diabwc=2,hp=2,rend=2,canc=2,msld=3,metacanc=6,aids=6)
 score<-as.numeric(as.matrix(z[names(weights)])%*%weights)
 if('charlson_noncancer'%in%names(co)){old<-co$charlson_noncancer[match(z$id,co$hospitalization_id)];stopifnot(all(is.na(old)|score==old))}
 z<-z |> dplyr::select(-canc,-metacanc);names(z)[names(z)=='id']<-'hospitalization_id'
 names(z)[names(z)!='hospitalization_id']<-paste0('charlson_',names(z)[names(z)!='hospitalization_id'])
 z
}
descriptive_enrich<-function(d,panel=NULL){
 unknown<-function(x)is.na(x)|!nzchar(trimws(x))|grepl('unknown|missing|refused|declined|not reported|unavailable|not collected',x,ignore.case=TRUE)
 if('sex_category'%in%names(d)){x<-tolower(trimws(d$sex_category));d$sex_group<-ifelse(unknown(x),'Unknown',ifelse(x%in%c('male','m'),'Male',ifelse(x%in%c('female','f'),'Female','Other recorded sex')))}
 if('race_category'%in%names(d)){x<-tolower(trimws(d$race_category));d$race_group<-ifelse(unknown(x),'Unknown',ifelse(grepl('black|african',x),'Black',ifelse(grepl('white',x),'White','Other known race')))}
 if('ethnicity_category'%in%names(d)){x<-tolower(trimws(d$ethnicity_category));d$ethnicity_group<-ifelse(unknown(x),'Unknown/unspecified',ifelse(grepl('not hispanic|non.?hispanic',x),'Not Hispanic/Latino',ifelse(grepl('hispanic|latino',x),'Hispanic/Latino','Unknown/unspecified')))}

 d$ward_to_icu_hours<-hours_between(d$t0,d$admission_dttm);d$ward_to_icu_hours[d$ward_to_icu_hours<0]<-NA_real_
 if('discharge_dttm'%in%names(d)){d$hospital_los_days<-hours_between(d$discharge_dttm,d$admission_dttm)/24;d$hospital_los_days[d$hospital_los_days<0]<-NA_real_}
 d$whole_stay_outcome<-ifelse(d$valid_outcome,d$event_type,'unknown')
 for(pair in list(c('whole_stay_death','death'),c('whole_stay_hospice','hospice'),c('whole_stay_composite','composite'))){v<-if(pair[2]=='composite')d$event_type%in%c('death','hospice') else d$event_type==pair[2];d[[pair[1]]]<-ifelse(d$valid_outcome,v,NA)}
 if('time_from_landmark'%in%names(d)){
  d$outcome28<-ifelse(d$time_from_landmark<=28,d$event_type,'no_terminal_event_by_day28')
  d$death28<-d$time_from_landmark<=28&d$event_type=='death';d$hospice28<-d$time_from_landmark<=28&d$event_type=='hospice';d$composite28<-d$death28|d$hospice28
 }
 for(pair in list(c('any_low_o2','hours_LOW_O2'),c('any_hfnc','hours_HFNC'),c('any_cpap','hours_CPAP'),c('any_nippv','hours_NIPPV')))if(pair[2]%in%names(d))d[[pair[1]]]<-d[[pair[2]]]>0
 if(!is.null(panel)){
  z<-panel |> dplyr::arrange(h) |> dplyr::group_by(hospitalization_id) |> dplyr::summarise(first_hour_support=dplyr::first(support),hypox_positive_hours=sum(hypox=='H1'),hyper_positive_hours=sum(hyper=='C1'),mixed_positive_hours=sum(hypox=='H1'&hyper=='C1'),hypox_unknown_hours=sum(hypox=='HU'),hyper_unknown_hours=sum(hyper=='CU'),.groups='drop')
  z$any_hypox<-z$hypox_positive_hours>0;z$any_hyper<-z$hyper_positive_hours>0;z$any_mixed<-z$mixed_positive_hours>0
  d<-dplyr::left_join(d,z,by='hospitalization_id')
 }
 d
}
descriptive_bins<-function(v){
 if(v=='age_years')return(list(edges=c(0,18,seq(20,120,5),Inf),integer=FALSE))
 if(v%in%c('charlson_noncancer','sofa_total','sofa_observed_domains','apache_adapted_total'))return(list(edges=c(-Inf,seq(-.5,99.5,1),Inf),integer=TRUE))
 if(grepl('^hours_|_hours$',v)&&v!='ward_to_icu_hours')return(list(edges=c(-Inf,seq(-.5,72.5,1),Inf),integer=TRUE))
 if(v=='ward_to_icu_hours')return(list(edges=c(-Inf,0,6,12,24,48,72,120,168,336,720,Inf),integer=FALSE))
 if(v=='hospital_los_days')return(list(edges=c(-Inf,0,.25,.5,seq(1,30,1),45,60,90,180,365,Inf),integer=FALSE))
 list(edges=c(-Inf,seq(0,1,.05),Inf),integer=FALSE)
}
descriptive_summary<-function(d,cohort,groups='Overall',group_variable=NULL,min_cell=10){
 registry<-descriptive_registry();con<-list();cat<-list();hist<-list();avail<-list();den<-list()
 for(group in groups){
  z<-if(group=='Overall')d else d[d[[group_variable]]==group & !is.na(d[[group_variable]]),,drop=FALSE]
  den[[group]]<-tibble::tibble(cohort=cohort,group=group,n=nrow(z),site_distinct_patients_n=dplyr::n_distinct(z$patient_id),patient_count_additive=FALSE)
  for(i in seq_len(nrow(registry))){
   v<-registry$variable[i];applicable<-!(cohort=='eligible_cohort'&&registry$primary_only[i]);present<-v%in%names(z)&&applicable
   status<-if(!applicable)'not_calculated_for_this_cohort' else if(!present)'not_supplied_or_not_derived' else if(!nrow(z))'structural_zero_group' else 'calculated'
   avail[[paste(group,v)]]<-tibble::tibble(cohort=cohort,group=group,variable=v,status=status)
   if(!present)next
   if(registry$type[i]=='continuous'){
    x<-num(z[[v]]);x[!is.finite(x)]<-NA_real_;obs<-x[!is.na(x)];n<-length(obs);missing<-sum(is.na(x));q<-if(n)unname(quantile(obs,c(.25,.5,.75))) else rep(NA_real_,3)
    con[[paste(group,v)]]<-tibble::tibble(cohort=cohort,group=group,variable=v,denominator_n=nrow(z),observed_n=n,missing_n=missing,sum=if(n)sum(obs) else NA_real_,sum_squares=if(n)sum(obs^2) else NA_real_,mean=if(n)mean(obs) else NA_real_,sd=if(n>1)sd(obs) else NA_real_,q1=q[1],median=q[2],q3=q[3],quantiles='site_only_do_not_average')
    b<-descriptive_bins(v);h<-graphics::hist(obs,breaks=b$edges,plot=FALSE,right=FALSE,include.lowest=TRUE)$counts
    hist[[paste(group,v)]]<-tibble::tibble(cohort=cohort,group=group,variable=v,bin=seq_along(h),lower=b$edges[-length(b$edges)],upper=b$edges[-1],n=as.numeric(h),integer_exact=b$integer,observed_n=length(obs),missing_n=missing)
   }else{
    value<-as.character(z[[v]]);miss<-is.na(value)|!nzchar(trimws(value))|grepl('^(unknown|missing)([/ _-]|$)|^(not reported|declined|refused|unavailable)$',tolower(trimws(value)))
    value[is.na(value)|!nzchar(trimws(value))]<-'MISSING'
    q<-tibble::tibble(category=value,category_missing=miss) |> dplyr::count(category,category_missing,name='n')
    q$cohort<-cohort;q$group<-group;q$variable<-v;q$denominator_n<-nrow(z);q$observed_n<-sum(!miss);q$missing_n<-sum(miss);q$percent_total<-if(nrow(z))100*q$n/nrow(z) else NA_real_;q$partition_complete<-TRUE
    cat[[paste(group,v)]]<-q
   }
  }
 }
 list(continuous=dplyr::bind_rows(con),categorical=dplyr::bind_rows(cat),histograms=dplyr::bind_rows(hist),availability=dplyr::bind_rows(avail),denominators=dplyr::bind_rows(den))
}
descriptive_protect<-function(tables,min_cell){
 # Suppress across overall and phenotype margins of the same variable when any
 # linked component is small, including complementary observed/missing counts.
 for(k in c('continuous','categorical','histograms')){
  x<-tables[[k]];count<-intersect(c('denominator_n','observed_n','missing_n','n'),names(x))
  tables[[k]]<-protect_aggregate(x,count,c('cohort','variable'),min_cell)
 }
 tables$denominators<-protect_aggregate(tables$denominators,c('n','site_distinct_patients_n'),'cohort',min_cell)
 tables
}
capture_descriptive_tables<-function(prepared,primary,cfg,root,panel=NULL,phenotypes=NULL,apache=NULL,consensus_id=NULL){
 if(dir.exists(file.path(root,'share/descriptives')))stop('Descriptor package exists; choose fresh output')
 co<-prepared$cohort;flags<-descriptive_comorbidities(co,prepared$tables$hospital_diagnosis)
 add<-function(d){d<-dplyr::left_join(d,flags,by='hospitalization_id');descriptive_enrich(d)}
 co<-add(co);d<-add(primary)
 if(!is.null(panel))d<-descriptive_enrich(d,panel)
 if(!is.null(apache))d<-dplyr::left_join(d,apache |> dplyr::select(hospitalization_id,apache_adapted_total),by='hospitalization_id')
 groupvar<-if(!is.null(phenotypes))'phenotype' else 'local_cluster'
 groups<-if(!is.null(phenotypes))c('Overall',phenotypes,'I') else c('Overall',as.character(sort(unique(d$local_cluster))))
 tables<-list(eligible=descriptive_summary(co,'eligible_cohort'),landmark=descriptive_summary(d,'landmark_cohort'),phenotypes=descriptive_summary(d,if(groupvar=='phenotype')'primary_phenotypes' else'provisional_local_clusters',groups,groupvar))
 all<-lapply(names(tables[[1]]),function(k)dplyr::bind_rows(lapply(tables,`[[`,k)));names(all)<-names(tables[[1]])
 protected<-descriptive_protect(all,cfg$analysis$min_cell)
 dir.create(file.path(root,'local/descriptives'),recursive=TRUE,showWarnings=FALSE);dest<-file.path(root,'share/descriptives');dir.create(dest,recursive=TRUE)
 for(k in names(all)){
  readr::write_csv(all[[k]],file.path(root,'local/descriptives',paste0(k,'_exact.csv')),na='')
  z<-protected[[k]];z$site<-cfg$site_name;z$run_id<-cfg$analysis$run_id;z$consensus_id<-if(is.null(consensus_id))'not_yet_shared' else consensus_id
  readr::write_csv(z,file.path(dest,paste0(k,'.csv')),na='')
 }
 dictionary<-descriptive_registry();dictionary$site<-cfg$site_name;readr::write_csv(dictionary,file.path(dest,'dictionary.csv'))
 # Verify local invariants before releasing a completion manifest.
 for(x in tables){
  z<-x$continuous;stopifnot(all(z$observed_n+z$missing_n==z$denominator_n))
  z<-x$categorical;if(nrow(z)){p<-z |> dplyr::group_by(cohort,group,variable) |> dplyr::summarise(sum_n=sum(n),den=dplyr::first(denominator_n),.groups='drop');stopifnot(all(p$sum_n==p$den))}
  z<-x$histograms;if(nrow(z)){p<-z |> dplyr::group_by(cohort,group,variable) |> dplyr::summarise(sum_n=sum(n),obs=dplyr::first(observed_n),.groups='drop');stopifnot(all(p$sum_n==p$obs))}
 }
 for(f in list.files(dest,pattern='csv$',full.names=TRUE)){z<-readr::read_csv(f,n_max=0,show_col_types=FALSE);stopifnot(!any(names(z)%in%c('patient_id','hospitalization_id','admission_dttm','discharge_dttm','t0','event_time','recorded_dttm')))}
 jsonlite::write_json(list(status='passed',site=cfg$site_name,run_id=cfg$analysis$run_id,consensus_id=consensus_id,grouping=groupvar,cohorts=c('eligible_cohort','landmark_cohort',if(groupvar=='phenotype')'primary_phenotypes' else'provisional_local_clusters'),continuous_pooling='observed n,sum,sum_squares; do not average site medians',histograms='fixed edges; exact quantiles only for integer bins; otherwise quantile bounds',categorical='complete category partitions with total/observed/missing denominators; unknown not zero',patients='site-distinct counts not additive across phenotypes or sites',outcomes='whole-stay descriptive and28d-after-landmark explicitly separate; LOS descriptive only',aggregate_policy=if(aggregate_exact())'exact aggregate counts; no disclosure suppression' else 'legacy linked suppression',source_id=prepared$source_id,source_md5=as.list(tools::md5sum('code/lib/descriptives.R')),files=as.list(setNames(unname(tools::md5sum(list.files(dest,full.names=TRUE))),list.files(dest)))),file.path(dest,'manifest.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
 invisible(dest)
}
descriptive_panel_from_sequence<-function(fit){
 sq<-fit$sequence;stopifnot(length(dim(sq))==2,nrow(sq)==length(fit$ids))
 states<-as.vector(sq);phys<-sub('^[^|]*[|]','',states)
 tibble::tibble(hospitalization_id=rep(fit$ids,ncol(sq)),h=rep(seq_len(ncol(sq))-1L,each=nrow(sq)),support=sub('[|].*$','',states),hypox=sub('_.*$','',phys),hyper=sub('^.*_','',phys))
}
