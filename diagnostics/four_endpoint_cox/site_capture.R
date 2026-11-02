# Add association endpoints to frozen assignments; patient records remain local.
source('diagnostics/federated/model_lib.R');source('diagnostics/federated/cox_lib.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=4)stop('Usage: site_capture.R SITE PRIMARY_ROOT REVISION_ROOT NEW_OUTPUT')
site<-a[1];primary<-a[2];revision<-a[3];out<-a[4]
if(dir.exists(out))stop('Choose a fresh output root')
variants<-c('primary','h24','h48','onset','support_only','low_flow','carry0','carry12','raw_support')
inputs<-list();audits<-list();shared<-file.path(out,'share');private<-file.path(out,'local')
for(variant in variants){
 root<-if(variant=='primary')primary else file.path(revision,'sensitivities',variant)
 home<-file.path(root,'sites',site);cache<-file.path(home,'local/full/outcomes.rds')
 old<-readRDS(cache);source_hash<-unname(tools::md5sum(cache));oldmeta<-readr::read_csv(file.path(home,'export/outcomes/full/model_metadata.csv'),show_col_types=FALSE)
 cm<-jsonlite::fromJSON(file.path(root,'coordinator/full/frozen_manifest.json'))
 phenotypes<-rownames(readRDS(file.path(root,'coordinator/full/prototypes.rds')))
 x<-fed_cox_associations(old$data,phenotypes,site,'full',cm$consensus_id)
 # All previously fitted endpoints must reproduce exactly, including robust SEs.
 maxdiff<-0
 for(key in names(old$cox)){
  before<-old$cox[[key]];after<-x$fits[[key]]
  if(inherits(before,'coxph')){
   stopifnot(inherits(after,'coxph'),before$n==after$n,before$nevent==after$nevent,identical(names(coef(before)),names(coef(after))))
   delta<-c(abs(coef(before)-coef(after)),abs(vcov(before)-vcov(after)));delta<-delta[is.finite(delta)]
   if(length(delta))maxdiff<-max(maxdiff,max(delta))
  }else stopifnot(inherits(after,'error'))
 }
 stopifnot(maxdiff<1e-9,identical(old$data$hospitalization_id,x$data$hospitalization_id),identical(old$data$phenotype,x$data$phenotype),unname(tools::md5sum(cache))==source_hash)
 for(endpoint in setdiff(cox_endpoints(),'live_discharge')){
  before<-oldmeta[oldmeta$endpoint==endpoint,];after<-x$metadata[x$metadata$endpoint==endpoint,]
  normalize_status<-function(s)trimws(ifelse(is.na(s),'',s))
  stopifnot(identical(normalize_status(before$status),normalize_status(after$status)),identical(as.numeric(before$events_lower),as.numeric(after$events_lower)))
 }
 x$metadata$status<-trimws(x$metadata$status);x$metadata$status[is.na(x$metadata$status)|!nzchar(x$metadata$status)]<-'model_nonestimable_nonfinite_coefficients'
 diagnosis<-supplement_model_diagnostics(list(data=x$data,cox=x$fits),x$metadata)
 stopifnot(!any(diagnosis$table$status%in%c('diagnostic_reconstruction_mismatch','diagnostic_refit_failed')))
 piece<-supplement_piecewise_associations(x$data)
 for(endpoint in cox_endpoints())stopifnot(sum(piece$metadata$events_n[piece$metadata$endpoint==endpoint])==sum(cox_event(x$data,endpoint)))
 partition<-x$data |> dplyr::group_by(phenotype) |> dplyr::summarise(n=dplyr::n(),death=sum(cox_event(dplyr::pick(dplyr::everything()),'death')),hospice=sum(cox_event(dplyr::pick(dplyr::everything()),'hospice')),live_discharge=sum(cox_event(dplyr::pick(dplyr::everything()),'live_discharge')),.groups='drop')
 partition$no_terminal_event_by28<-partition$n-partition$death-partition$hospice-partition$live_discharge
 stopifnot(all(partition$no_terminal_event_by28>=0))
 tables<-list(associations=x$associations,model_metadata=x$metadata,association_covariances=x$covariances,event_partitions=partition,model_PH_and_denominators=diagnosis$table,model_PH_terms=diagnosis$term_tests,exploratory_piecewise_associations=piece$associations,exploratory_piecewise_metadata=piece$metadata)
 for(name in names(tables)){
  z<-tables[[name]];z$site<-rep(site,nrow(z));z$analysis_variant<-rep(variant,nrow(z));z$consensus_id<-rep(cm$consensus_id,nrow(z))
  fed_write(z,file.path(shared,variant,paste0(name,'.csv')))
 }
 dir.create(file.path(private,variant),recursive=TRUE,showWarnings=FALSE)
 saveRDS(list(data=x$data,cox=x$fits,diagnostics=diagnosis,piecewise=piece),file.path(private,variant,'models.rds'))
 audits[[variant]]<-tibble::tibble(site=site,analysis_variant=variant,n=nrow(x$data),original_endpoints_max_abs_difference=maxdiff,assignments_unchanged=TRUE,source_cache_unchanged=TRUE,models=nrow(x$metadata),PH_reconstruction_passed=TRUE)
 inputs[[variant]]<-list(consensus_id=cm$consensus_id,outcomes_md5=source_hash)
 cat(site,variant,'four endpoints complete; original effects reproduced\n');flush.console()
}
fed_write(dplyr::bind_rows(audits),file.path(shared,'audit.csv'))
f<-list.files(shared,recursive=TRUE,full.names=FALSE)
fed_json(list(status='passed',site=site,primary_consensus_id=inputs$primary$consensus_id,variants=variants,endpoints=cox_endpoints(),horizon_days=28,model_specification=model_spec(),primary_stratum='I retained with separate baseline hazard; P1 reference',live_discharge='Live non-hospice discharge; death and hospice censored at competing exit',original_endpoints_reproduced=TRUE,assignments_unchanged=TRUE,models_refit='association models only; frozen prediction models unchanged',inputs=inputs,source_md5=as.list(tools::md5sum(c('diagnostics/four_endpoint_cox/site_capture.R','diagnostics/federated/cox_lib.R','code/lib/models.R','code/lib/site_supplement.R'))),files=as.list(setNames(unname(tools::md5sum(file.path(shared,f))),f))),file.path(shared,'manifest.json'))
