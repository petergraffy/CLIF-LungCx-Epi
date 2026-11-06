# Data derivation and fixed-template classification only. No clustering is run.
source('code/lib/load.R')
s2_sequence<-function(panel,field){
 wide<-panel |> dplyr::select(hospitalization_id,h,state=dplyr::all_of(field)) |> tidyr::pivot_wider(names_from=h,values_from=state,names_prefix='h') |> dplyr::arrange(hospitalization_id)
 cols<-paste0('h',sort(unique(panel$h)))
 list(ids=wide$hospitalization_id,sequence=as.matrix(wide[cols]),field=field,method='fixed_template_assignment_only')
}
s2_raw_inputs<-function(cfg){
 prepared<-prepare_revision(cfg);root<-prepared$root;all<-prepared$cohort;tables<-prepared$tables;spec<-analysis_specifications(cfg)
 write_share(spec,root,'analysis_plan',cfg)
 detection<-first_detected_arf(all,tables,cfg);saveRDS(detection,file.path(root,'local/arf_detection.rds'))
 write_share(detection |> dplyr::count(arf_detection_status,arf_initial_evidence,name='n'),root,'arf_detection_audit',cfg)
 analyses<-list()
 for(i in seq_len(nrow(spec))){
  window<-spec$window[i];variant<-spec$variant[i];id<-spec$analysis_id[i]
  d<-if(variant=='arf_onset')arf_onset_cohort(all,detection) else landmark_cohort(all,window)
  write_share(d |> dplyr::count(landmark_reason,name='n') |> dplyr::mutate(analysis_id=id),root,paste0(id,'_landmark_flow'),cfg)
  d<-if(variant=='absorbing')d[d$valid_outcome,] else d[d$landmark_reason=='eligible',]
  if(!nrow(d))stop('No eligible local observations for ',id,'. No clustering or template changes were attempted.')
  panel<-build_panel(d,tables,window,cfg,absorbing=variant=='absorbing',estimate_low_flow=variant=='low_flow_fio2',carry=switch(variant,arf_carry_0=0,arf_carry_12=12,NULL))
  d<-d |> dplyr::left_join(support_summaries(panel),by='hospitalization_id') |> dplyr::left_join(code_status_features(d,tables$code_status,window,cfg),by='hospitalization_id') |> dplyr::left_join(vasopressor_features(d,tables$medication_admin_continuous,window),by='hospitalization_id')
  d$sequence_reason<-'included';field<-if(variant=='support_only')'support' else 'joint'
  analyses[[id]]<-list(data=d,panel=panel,fit=s2_sequence(panel,field),window=window,variant=variant,specification=spec[i,])
 }
 saveRDS(list(analyses=analyses,parameters=cfg$analysis,source_id=code_fingerprint(),clustering_performed=FALSE),file.path(root,'local/derivation.rds'))
 write_share(tibble::tibble(analysis_id=names(analyses),status='data_rebuilt_no_clustering'),root,'derivation_status',cfg)
 invisible(analyses)
}
s2_waterfall_inputs<-function(cfg,out){
 im<-jsonlite::fromJSON(file.path(out,'local/input_manifest.json'));raw<-readRDS(file.path(im$baseline_root,'local/derivation.rds'));prepared<-readRDS(file.path(im$baseline_root,'local/prepared.rds'));results<-list()
 dir.create(file.path(out,'aggregate'),recursive=TRUE,showWarnings=FALSE)
 for(id in names(raw$analyses)){
  a<-raw$analyses[[id]];job<-if(id%in%c('h72_absorbing','h24_landmark_sensitivity','h48_landmark_sensitivity','h72_arf_onset'))id else 'h72_primary'
  rs<-as.data.frame(arrow::read_parquet(file.path(out,'local',paste0(job,'_waterfall.parquet'))));rs$recorded_dttm<-utc(rs$recorded_dttm);rs<-rs[rs$hospitalization_id%in%a$data$hospitalization_id,]
  limits<-a$data |> dplyr::transmute(hospitalization_id,start=pmax(admission_dttm,t0-cfg$analysis$support_carry_hours*3600),end=pmin(t0+a$window*3600,event_time))
  bounded<-rs |> dplyr::left_join(limits,by='hospitalization_id');stopifnot(all(bounded$recorded_dttm>=bounded$start & bounded$recorded_dttm<bounded$end))
  tables<-prepared$tables;tables$respiratory_support<-rs
  panel<-build_panel(a$data,tables,a$window,cfg,absorbing=a$variant=='absorbing',estimate_low_flow=a$variant=='low_flow_fio2',carry=switch(a$variant,arf_carry_0=0,arf_carry_12=12,NULL))
  old<-a$panel |> dplyr::transmute(hospitalization_id,h,old_support=support,old_physiology=physiology,old_joint=joint)
  panel<-panel |> dplyr::left_join(old,by=c('hospitalization_id','h'))
  direct<-rs |> dplyr::filter(raw_device_observed_at_timestamp) |> dplyr::left_join(a$data |> dplyr::select(hospitalization_id,t0),by='hospitalization_id') |> dplyr::mutate(h=floor(hours_between(recorded_dttm,t0)),support=support_state(device_category))
  panel$support_observed<-paste(panel$hospitalization_id,panel$h,panel$support)%in%paste(direct$hospitalization_id,direct$h,direct$support)&!panel$terminal
  panel$support_carried<-panel$support!='UNKNOWN'&!panel$support_observed&!panel$terminal
  stopifnot(nrow(panel)==a$window*nrow(a$data),!anyNA(panel$old_support))
  f<-support_summaries(panel);d<-a$data |> dplyr::select(-dplyr::any_of(setdiff(names(f),'hospitalization_id'))) |> dplyr::left_join(f,by='hospitalization_id')
  results[[id]]<-modifyList(a,list(data=d,panel=panel,fit=s2_sequence(panel,if(a$variant=='support_only')'support' else 'joint')))
 }
 a<-results$h72_primary;a$panel$physiology<-a$panel$old_physiology;a$panel$hypox<-sub('_.*','',a$panel$physiology);a$panel$hyper<-sub('.*_','',a$panel$physiology);a$panel$joint<-paste(a$panel$support,a$panel$physiology,sep='|')
 f<-support_summaries(a$panel);a$data<-a$data |> dplyr::select(-dplyr::any_of(setdiff(names(f),'hospitalization_id'))) |> dplyr::left_join(f,by='hospitalization_id');a$fit<-s2_sequence(a$panel,'joint');results$h72_waterfall_support_with_raw_arf<-a
 saveRDS(list(analyses=results,clustering_performed=FALSE),file.path(out,'local/waterfall_derivation.rds'))
}
