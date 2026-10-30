source('diagnostics/revision20261022/score_lib.R')
site<-commandArgs(trailingOnly=TRUE)[1];root<-rev_root();home<-file.path(root,'sites',site);checks<-list();abl<-list()
original<-readRDS(file.path(old_root(),'sites',site,'local/full/outcomes.rds'))$data
for(v in c('h24','h48','onset','support_only','low_flow','carry0','carry12','raw_support')){
 Sys.setenv(CLIF_REVISION_VARIANT=v);source('diagnostics/revision20261022/variant_contract.R')
 a<-readRDS(file.path(root,'sensitivities',v,'sites',site,'local/discovery.rds'));d<-a$data
 stopifnot(nrow(a$panel)==nrow(d)*REV_WINDOW,all(d$time_from_landmark>0),!any(a$panel$terminal))
 modal<-fed_modal(a$panel);stopifnot(setequal(unique(modal$bin[modal$domain=='support']),0:11),all(modal$state[modal$domain=='joint']%in%fed_domains()$joint))
 for(version in c('full','loo_UCMC','loo_NU','loo_RUSH')){
 z<-readRDS(file.path(root,'sensitivities',v,'sites',site,'local',version,'outcomes.rds'))$data
 stopifnot(nrow(z)==nrow(d),!anyNA(z$phenotype),setequal(z$hospitalization_id,d$hospitalization_id),all(z$followup>0 & z$followup<=28),all(z$phenotype[z$hours_UNKNOWN==REV_WINDOW]=='I'))
 }
 if(v%in%c('support_only','raw_support')){
 z<-readRDS(file.path(root,'sensitivities',v,'sites',site,'local/full/outcomes.rds'))$data;stopifnot(setequal(original$hospitalization_id,z$hospitalization_id));z<-z[match(original$hospitalization_id,z$hospitalization_id),]
 ar<-adjusted_rand(as.integer(factor(original$phenotype)),as.integer(factor(z$phenotype)))
 ix<-original$phenotype!='I' & z$phenotype!='I';ar2<-adjusted_rand(as.integer(factor(original$phenotype[ix])),as.integer(factor(z$phenotype[ix])))
 abl[[v]]<-tibble::tibble(site=site,analysis=v,adjusted_rand_all=ar,adjusted_rand_both_classified=ar2,interpretation='Agreement of partitions, not label accuracy; independently learned shared templates')
 }
 checks[[v]]<-list(n=nrow(d),window=REV_WINDOW,all_retained=TRUE,all_twelve_bins_present=TRUE,no_terminal_or_prelandmark_outcomes=TRUE)
}
# Independently reproduce score predictions using only the frozen other-site coefficients.
source('diagnostics/revision20261022/score_lib.R')
version<-paste0('loo_',site);m<-jsonlite::fromJSON(file.path(root,'coordinator',version,'score_predictor.json'));stopifnot(!site%in%m$training_sites,!any(grepl(paste0('/scores/',site,'/'),m$input_files)))
b<-readr::read_csv(file.path(root,'coordinator',version,'score_coefficients.csv'),show_col_types=FALSE)
d<-readRDS(file.path(old_root(),'sites',site,'local',version,'outcomes.rds'))$data
pv<-readRDS(file.path(home,'local/score_transport.rds'))
for(cohort in score_spec()$cohorts)for(endpoint in score_spec()$endpoints){
 z<-if(cohort=='measured_sofa')d[!is.na(d$sofa_total),] else d;pr<-pv[[paste(cohort,endpoint)]];stopifnot(identical(z$hospitalization_id,pr$ids))
 for(kind in score_spec()$models){xx<-score_matrix(z,kind,cohort,revision_phenotypes(version));bb<-b[b$cohort==cohort & b$endpoint==endpoint & b$model==kind,];beta<-setNames(bb$estimate,bb$term);pred<-plogis(drop(xx%*%beta[colnames(xx)]));stopifnot(max(abs(pred-pr$predictions[[kind]]))<1e-12)}
}
# Non-overlapping original fixed artifacts remain unchanged.
ledger<-jsonlite::fromJSON(file.path(old_root(),'reports/frozen_artifact_ledger.json'),simplifyVector=FALSE)$ledger
for(a in ledger)stopifnot(unname(tools::md5sum(a$path))==a$md5)
fed_write(dplyr::bind_rows(abl),file.path(home,'export/partition_agreement.csv'))
fed_json(list(site=site,status='passed',sensitivities=checks,score_comparison_same_patients=TRUE,independent_frozen_prediction_check=TRUE,original_primary_artifacts_unchanged=TRUE,all_primary_cases_retained_in_missing_domain_sensitivity=TRUE),file.path(home,'export/revision_audit.json'))
cat(site,' revision local audits passed\n')
