source('diagnostics/federated/comparator_lib.R')
rev_root<-function()normalizePath(Sys.getenv('CLIF_REVISION_OUTPUT',unset='output/run_revision_gaps_20261022'),mustWork=FALSE)
old_root<-function()normalizePath('output/run_federated_rehearsal_20261022')
rev_write<-function(x,path)fed_write(x,file.path(rev_root(),path))
# Release interval-valued categorical counts; linked package still needs review.
band_table<-function(d,groups,variable,value){
 d$category<-as.character(value);d$category[is.na(d$category)]<-'unknown'
 z<-d |> dplyr::count(dplyr::across(dplyr::all_of(groups)),category,name='count')
 z$variable<-variable;z$n_lower<-10*floor(z$count/10);z$n_upper<-z$n_lower+9;z$count<-NULL;z
}
revision_input<-function(site){
 variant<-Sys.getenv('CLIF_REVISION_VARIANT');wf<-readRDS(paste0('output/run_',site,'_waterfall_20261022/local/waterfall_derivation.rds'))$analyses
 m<-jsonlite::fromJSON(paste0('output/run_',site,'_waterfall_20261022/local/input_manifest.json'))
 raw<-readRDS(file.path(m$baseline_root,'local/derivation.rds'))$analyses
 key<-switch(variant,h24='h24_landmark_sensitivity',h48='h48_landmark_sensitivity',onset='h72_arf_onset',support_only='h72_primary',low_flow='h72_low_flow_fio2',carry0='h72_arf_carry_0',carry12='h72_arf_carry_12',raw_support='h72_primary',stop('Unknown sensitivity'))
 a<-wf[[key]];sq<-raw[[key]]$fit$sequence
 if(raw[[key]]$fit$field=='support')stop('Need joint raw physiology')
 phys<-sub('^[^|]*\\|','',as.vector(sq));orig<-tibble::tibble(hospitalization_id=rep(raw[[key]]$fit$ids,ncol(sq)),h=rep(0:(ncol(sq)-1),each=nrow(sq)),raw_phys=phys)
 if(variant=='raw_support')a$panel$support<-a$panel$old_support
 a$panel<-a$panel |> dplyr::left_join(orig,by=c('hospitalization_id','h'))
 stopifnot(!anyNA(a$panel$raw_phys));a$panel$physiology<-a$panel$raw_phys
 a$panel$hypox<-sub('_.*','',a$panel$physiology);a$panel$hyper<-sub('.*_','',a$panel$physiology)
 a$panel$joint<-paste(a$panel$support,a$panel$physiology,sep='|')
 features<-support_summaries(a$panel)
 a$data<-a$data |> dplyr::select(-dplyr::any_of(c('phenotype','assignment_status','local_cluster',setdiff(names(features),'hospitalization_id')))) |> dplyr::left_join(features,by='hospitalization_id')
 if(variant=='support_only')a$panel$joint<-a$panel$support
 a
}
