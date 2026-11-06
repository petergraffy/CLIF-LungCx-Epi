# Collect a mixture of native and verified Python-fallback Stage 2 packets.
source('python_fallback/environment.R');clif_fallback_check_environment();source('scripts/deployment_lib.R');source('code/lib/central.R');source('python_fallback/verify_return.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)<3L)stop('Usage: central_collect.R BUNDLE NEW_OUTPUT SITE_RETURN [...]')
bundle<-a[1];out<-a[2];returns<-a[-c(1,2)];bm<-deploy_validate_inventory(bundle,'bundle_manifest.json');ms<-lapply(returns,deploy_validate_inventory)
sites<-vapply(ms,`[[`,character(1),'site');stopifnot(!anyDuplicated(sites),setequal(sites,unlist(bm$sites)))
bhash<-unname(tools::md5sum(file.path(bundle,'bundle_manifest.json')))
for(i in seq_along(ms)){m<-ms[[i]];verify_stage2_return(m,bm,bhash);deploy_check_csvs(returns[i])}
if(dir.exists(out))stop('Choose fresh collection output')
dir.create(out,recursive=TRUE)
# Retain complete site tables with explicit site labels for every supplement/figure.
allfiles<-unique(unlist(lapply(returns,function(p)list.files(p,pattern='[.]csv$',recursive=TRUE))))
for(rel in allfiles){z<-list();for(i in seq_along(returns)){f<-file.path(returns[i],rel);if(file.exists(f)){x<-readr::read_csv(f,show_col_types=FALSE);if(!'site'%in%names(x))x$site<-rep(sites[i],nrow(x));z[[sites[i]]]<-x}}
 combined<-dplyr::bind_rows(z);dest<-file.path(out,'site_tables',rel);dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE);readr::write_csv(combined,dest,na='')
}
# Same pooling estimator as the established association pipeline, without freezing predictors.
for(v in setdiff(DEPLOY_VARIANTS,'absorbing')){
 z<-dplyr::bind_rows(lapply(returns,function(p)readr::read_csv(file.path(p,'analyses',v,'outcomes/full/associations.csv'),show_col_types=FALSE)))
 pooled<-z |> dplyr::group_by(model,endpoint,phenotype,reference,effect_type) |> dplyr::group_modify(~pool_effect(.x)) |> dplyr::ungroup()
 f<-file.path(out,'pooled_associations',paste0(v,'.csv'));dir.create(dirname(f),recursive=TRUE,showWarnings=FALSE);readr::write_csv(pooled,f,na='')
}
for(name in c('sensitivity_associations','arf_associations')){
 z<-dplyr::bind_rows(lapply(returns,function(p)readr::read_csv(file.path(p,'review',paste0(name,'.csv')),show_col_types=FALSE)))
 groups<-intersect(c('subset','model','endpoint','phenotype','reference','term'),names(z))
 pooled<-z |> dplyr::group_by(dplyr::across(dplyr::all_of(groups))) |> dplyr::group_modify(~pool_effect(.x)) |> dplyr::ungroup()
 readr::write_csv(pooled,file.path(out,'pooled_associations',paste0(name,'.csv')),na='')
}
deploy_json(list(status='passed',sites=sites,primary_consensus_id=bm$primary_consensus_id,workflow=bm$workflow,patient_level_inputs_read=FALSE,parquet_backends=setNames(lapply(ms,function(m)if(is.null(m$python_fallback))'R_arrow' else 'pyarrow'),sites),prediction_validation='local patient-grouped cross-validation; pooled out-of-fold risk histograms; not external-site predictor validation'),file.path(out,'collection_manifest.json'))
