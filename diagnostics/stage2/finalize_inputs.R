source('scripts/deployment_lib.R')
a<-commandArgs(trailingOnly=TRUE);site<-a[1];root<-a[2];cfg<-deploy_config(Sys.getenv('CLIF_RUN_CONFIG'));wf<-file.path(root,'waterfall')
raw<-readRDS(file.path(cfg$work_root,'local/derivation.rds'));analyses<-list()
keys<-c(primary='h72_primary',h24='h24_landmark_sensitivity',h48='h48_landmark_sensitivity',onset='h72_arf_onset',support_only='h72_support_only',low_flow='h72_low_flow_fio2',carry0='h72_arf_carry_0',carry12='h72_arf_carry_12',absorbing='h72_absorbing')
for(v in names(keys)){
 c<-readRDS(file.path(deploy_analysis_root(root,v),'sites',site,'local/discovery.rds'));a<-raw$analyses[[keys[[v]]]]
 a$data<-c$data;a$panel<-c$panel;a$fit<-c$fit;a$consensus_id<-c$consensus_id;a$grouping<-'frozen_template';analyses[[keys[[v]]]]<-a
}
saveRDS(list(analyses=analyses,clustering_performed=FALSE),file.path(cfg$work_root,'local/stage2_derivation.rds'))
w<-readRDS(file.path(wf,'local/waterfall_derivation.rds'));w$analyses$h72_absorbing<-analyses$h72_absorbing;saveRDS(w,file.path(wf,'local/waterfall_derivation.rds'))
full<-analyses$h72_primary$data;bm<-deploy_read(file.path(commandArgs(trailingOnly=TRUE)[3],'bundle_manifest.json'));reference<-bm$stage1_reference_counts[[site]]
now<-c(eligible=nrow(readRDS(file.path(cfg$work_root,'local/prepared.rds'))$cohort),landmark=nrow(full))
z<-tibble::tibble(site=site,cohort=names(now),stage2_n=unname(now),stage1_reference_n=vapply(names(now),function(k)if(is.null(reference[[k]]))NA_real_ else reference[[k]],numeric(1)))
z$difference<-z$stage2_n-z$stage1_reference_n;z$status<-ifelse(is.na(z$difference),'reference_unavailable',ifelse(z$difference==0,'count_matches_stage1','count_differs_review_source_or_derivation'))
dir.create(file.path(root,'input_audit'),showWarnings=FALSE);readr::write_csv(z,file.path(root,'input_audit/cohort_rebuild_comparison.csv'))
deploy_json(list(site=site,mode='rebuild_from_clif',stage1_files_required=FALSE,clustering_performed=FALSE,frozen_phenotypes_relearned=FALSE,primary_consensus_id=bm$primary_consensus_id,timestamp_policy=cfg$source_timezone%||%'UTC',count_comparison='Counts only; matching counts do not establish identical patients.'),file.path(root,'input_audit/rebuild_manifest.json'))
