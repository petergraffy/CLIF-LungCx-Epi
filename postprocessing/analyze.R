source("scripts/deployment_lib.R")
source("code/lib/central.R")
source("python_fallback/verify_return.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Usage: analyze.R BUNDLE EXISTING_OUTPUT SITE_RETURN [...]")
bundle <- args[[1L]]
out <- args[[2L]]
paths <- args[-c(1L, 2L)]
bm <- deploy_validate_inventory(bundle, "bundle_manifest.json")
bhash <- unname(tools::md5sum(file.path(bundle, "bundle_manifest.json")))
release <- deploy_read("python_fallback/release.json")
manifests <- lapply(paths, deploy_validate_inventory)
sites <- vapply(manifests, `[[`, character(1), "site")
stopifnot(!anyDuplicated(sites), all(sites %in% unlist(bm$sites)))
read <- function(path, rel) readr::read_csv(file.path(path, rel), show_col_types = FALSE)

rows<-list();metadata<-assignment<-agreement<-continuous<-categorical<-cohorts<-localauc<-assoc<-list()
for (i in seq_along(paths)) {
 path <- paths[[i]]
 m <- manifests[[i]]
 site <- sites[[i]]
 deploy_check_csvs(path)
 cat('Auditing ',site,'\n',sep='');flush.console()
 stopifnot(m$status=='complete',m$stage=='stage2',m$workflow==bm$workflow,
           m$r_lock_md5==bm$r_lock_md5, m$primary_consensus_id==bm$primary_consensus_id)
 issued<-tryCatch({verify_stage2_return(m,bm,bhash);m$r_version=='4.4.2'},error=function(e)FALSE)
 fallback<-!is.null(m$python_fallback)
 expected<-if(fallback) release$effective_source_files_md5 else release$parent_source_files_md5
 integrity <- 'verified'
 changed<-names(expected)[vapply(names(expected),function(k)!identical(m$source_files_md5[[k]],expected[[k]]),logical(1))]
 # Keep the already reviewed UCSF operational difference explicit in the audit.
 ucsf_operational_difference <- site == "UCSF" && !fallback &&
   identical(m$parameter_bundle_md5, bhash) && identical(m$r_version, "4.4.3") &&
   setequal(names(expected), names(m$source_files_md5)) &&
   setequal(changed, c("scripts/deployment_lib.R", "scripts/setup.R", "scripts/site.R"))
 if (!issued && !ucsf_operational_difference) stop("Unrecognized Stage 2 contract for ", site)
 if (!issued) warning("Including the provisional UCSF return; see packet_audit.csv", call. = FALSE)
 rows[[site]]<-data.frame(site=site,files=length(m$files),inventory=integrity,same_bundle=issued || identical(m$parameter_bundle_md5,bhash),issued_contract=issued,backend=if(fallback)'pyarrow' else 'R Arrow',same_consensus=identical(m$primary_consensus_id,bm$primary_consensus_id),source_exact=identical(m$source_hash,bm$source_hash),r_version=m$r_version,changed_files=paste(changed,collapse=';'))
 metadata[[site]]<-read(path,'analyses/primary/outcomes/full/model_metadata.csv')
 assignment[[site]]<-read(path,'analyses/primary/outcomes/full/assignment_counts.csv')
 agreement[[site]]<-read(path,'profiles/phenotype_transport_agreement.csv')
 continuous[[site]]<-read(path,'descriptives/continuous.csv')
 categorical[[site]]<-read(path,'descriptives/categorical.csv')
 cohorts[[site]]<-read(path,'input_audit/cohort_rebuild_comparison.csv')
 localauc[[site]]<-read(path,'prediction/local_validation.csv')
 assoc[[site]]<-read(path,'analyses/primary/outcomes/full/associations.csv')
}
audit<-dplyr::bind_rows(rows);readr::write_csv(audit,file.path(out,'packet_audit.csv'))
stopifnot(all(audit$inventory=='verified'),all(audit$same_bundle),all(audit$same_consensus))
for(k in c('metadata','assignment','agreement','continuous','categorical','cohorts','localauc','assoc'))assign(k,dplyr::bind_rows(get(k)))
for(k in c('metadata','assignment','agreement','continuous','categorical','cohorts','localauc','assoc'))readr::write_csv(get(k),file.path(out,paste0('combined_',k,'.csv')))
site_counts<-cohorts |> dplyr::select(site,cohort,stage2_n) |> tidyr::pivot_wider(names_from=cohort,values_from=stage2_n)
pheno<-assignment |> dplyr::group_by(site,phenotype) |> dplyr::summarise(n=sum(n),.groups='drop') |> tidyr::pivot_wider(names_from=phenotype,values_from=n,values_fill=0)
site_counts<-dplyr::left_join(site_counts,pheno,by='site')
agree<-agreement |> dplyr::group_by(site) |> dplyr::summarise(n_total=sum(n),agreed=sum(n[full==loo_aligned]),percent_agreed=100*agreed/n_total,.groups='drop')
site_counts<-dplyr::left_join(site_counts,agree |> dplyr::select(site,percent_agreed),by='site')
readr::write_csv(site_counts,file.path(out,'site_counts.csv'))
# Exact day28 event counts from the descriptor partitions, including information-insufficient stays.
event<-categorical |> dplyr::filter(cohort=='primary_phenotypes',variable=='outcome28') |> dplyr::group_by(group,category) |> dplyr::summarise(n=sum(n),.groups='drop') |> tidyr::pivot_wider(names_from=category,values_from=n,values_fill=0)
readr::write_csv(event,file.path(out,'phenotype_day28_counts.csv'))
# REML modified Hartung-Knapp, as in the study's existing pooling code.
pooled<-assoc |> dplyr::group_by(model,endpoint,phenotype,reference,effect_type) |> dplyr::group_modify(~pool_effect(.x)) |> dplyr::ungroup()
readr::write_csv(pooled,file.path(out,'pooled_primary.csv'))
exact_sites<-audit$site[audit$issued_contract]
strict<-assoc |> dplyr::filter(site%in%exact_sites) |> dplyr::group_by(model,endpoint,phenotype,reference,effect_type) |> dplyr::group_modify(~pool_effect(.x)) |> dplyr::ungroup()
readr::write_csv(strict,file.path(out,'pooled_primary_issued_contract_sites.csv'))
score<-continuous |> dplyr::filter(cohort=='landmark_cohort',group=='Overall',variable%in%c('sofa_total','apache_adapted_total')) |> dplyr::select(site,variable,denominator_n,observed_n,missing_n)
readr::write_csv(score,file.path(out,'severity_completeness.csv'))
flags<-metadata |> dplyr::filter(status!='ok');readr::write_csv(flags,file.path(out,'primary_model_flags.csv'))
cat('SITE COUNTS\n');print(site_counts,width=150)
cat('DAY28 OUTCOMES\n');print(event,width=150)
cat('CORE COX\n');print(pooled |> dplyr::filter(model=='core') |> dplyr::select(endpoint,phenotype,k,hazard_ratio,hr_low,hr_high,i2_percent),n=20,width=150)
cat('SEVERITY COMPLETENESS\n');print(score,n=20)
cat('MODEL FLAGS\n');print(flags |> dplyr::select(site,model,endpoint,status),n=30,width=160)

for(v in setdiff(DEPLOY_VARIANTS,'absorbing')){
  z<-dplyr::bind_rows(lapply(paths,read,rel=file.path('analyses',v,'outcomes/full/associations.csv')))
  fit<-z |> dplyr::group_by(model,endpoint,phenotype,reference,effect_type) |> dplyr::group_modify(~pool_effect(.x)) |> dplyr::ungroup()
  dir.create(file.path(out,'pooled_associations'),showWarnings=FALSE)
  readr::write_csv(fit,file.path(out,'pooled_associations',paste0(v,'.csv')))
}
for(name in c('sensitivity_associations','arf_associations')){
  z<-dplyr::bind_rows(lapply(paths,read,rel=file.path('review',paste0(name,'.csv'))))
  groups<-intersect(c('subset','model','endpoint','phenotype','reference','term'),names(z))
  fit<-z |> dplyr::group_by(dplyr::across(dplyr::all_of(groups))) |> dplyr::group_modify(~pool_effect(.x)) |> dplyr::ungroup()
  readr::write_csv(fit,file.path(out,'pooled_associations',paste0(name,'.csv')))
}

# Pool the exploratory interval models with the same estimator as the main Cox models.
piecewise <- dplyr::bind_rows(lapply(paths, read,
  rel = "four_endpoint_cox/primary/exploratory_piecewise_associations.csv"))
piecewise <- piecewise |>
  dplyr::filter(status == "ok") |>
  dplyr::group_by(interval_start_day, interval_end_day, endpoint, phenotype, reference) |>
  dplyr::group_modify(~ pool_effect(.x)) |>
  dplyr::ungroup()
readr::write_csv(piecewise, file.path(out, "pooled_associations", "piecewise.csv"))
