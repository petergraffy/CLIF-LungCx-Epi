# Aggregate-only Table1/Table2 assembly; no site clinical/cache access.
source('code/lib/load.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)<2)stop('Usage: collect_descriptives.R NEW_OUTPUT SITE_RETURN [...]')
out<-a[1];returns<-a[-1];if(dir.exists(out))stop('Choose fresh descriptor collection')
paths<-file.path(returns,'descriptives');m<-lapply(paths,function(p)jsonlite::fromJSON(file.path(p,'manifest.json')))
sites<-vapply(m,function(x)x$site,character(1));stopifnot(!anyDuplicated(sites),all(vapply(m,function(x)x$grouping=='phenotype',logical(1))),length(unique(vapply(m,function(x)x$consensus_id,character(1))))==1)
for(i in seq_along(paths)){
 hashes<-unlist(m[[i]]$files,use.names=TRUE);stopifnot(all(unname(tools::md5sum(file.path(paths[i],names(hashes))))==unname(hashes)))
}
dir.create(out,recursive=TRUE);all<-list()
for(k in c('continuous','categorical','histograms','availability','denominators','dictionary')){
 z<-dplyr::bind_rows(lapply(seq_along(paths),function(i){
  x<-readr::read_csv(file.path(paths[i],paste0(k,'.csv')),show_col_types=FALSE)
  stopifnot(all(x$site==sites[i]))
  if('consensus_id'%in%names(x))stopifnot(all(x$consensus_id==m[[i]]$consensus_id))
  x
 }))
 stopifnot(!any(names(z)%in%c('patient_id','hospitalization_id','t0','event_time','recorded_dttm','admission_dttm','discharge_dttm')))
 all[[k]]<-z;readr::write_csv(z,file.path(out,paste0(k,'_by_site.csv')),na='')
}
# Zero-case groups can contribute zero only when explicitly declared structural.
expected<-function(cohort,group,variable){
 z<-all$availability |> dplyr::filter(.data$cohort==.env$cohort,.data$group==.env$group,.data$variable==.env$variable)
 valid<-z$site[z$status%in%c('calculated','structural_zero_group')]
 length(unique(valid))==length(sites)
}
keys<-c('cohort','group','variable');continuous<-list()
for(k in unique(do.call(paste,c(all$continuous[keys],sep='|')))){
 z<-all$continuous[do.call(paste,c(all$continuous[keys],sep='|'))==k,];first<-z[1,]
 safe<-expected(first$cohort,first$group,first$variable)&&length(unique(z$site))==length(sites)&&all(z$cell_status=='observed')
 n<-if(safe)sum(z$observed_n) else NA_real_;miss<-if(safe)sum(z$missing_n) else NA_real_;den<-if(safe)sum(z$denominator_n) else NA_real_
 s<-if(safe)sum(ifelse(z$observed_n==0,0,z$sum)) else NA_real_;ss<-if(safe)sum(ifelse(z$observed_n==0,0,z$sum_squares)) else NA_real_
 continuous[[k]]<-tibble::tibble(cohort=first$cohort,group=first$group,variable=first$variable,denominator_n=den,observed_n=n,missing_n=miss,mean=if(safe&&n>0)s/n else NA_real_,sd=if(safe&&n>1)sqrt(max(0,(ss-s*s/n)/(n-1))) else NA_real_,status=if(!safe)'incomplete_or_suppressed_site_partition' else if(n==0)'all_missing' else'observed_all_sites')
}
pool_cont<-dplyr::bind_rows(continuous);readr::write_csv(pool_cont,file.path(out,'pooled_continuous.csv'),na='')
categorical<-list()
for(k in unique(do.call(paste,c(all$categorical[keys],sep='|')))){
 z<-all$categorical[do.call(paste,c(all$categorical[keys],sep='|'))==k,];first<-z[1,]
 safe<-expected(first$cohort,first$group,first$variable)&&all(z$cell_status=='observed')
 # A site may have no rows only for a declared zero-case group.
 for(site in setdiff(sites,unique(z$site))){q<-all$denominators |> dplyr::filter(.data$cohort==first$cohort,.data$group==first$group,.data$site==.env$site);safe<-safe&&nrow(q)==1&&!is.na(q$n)&&q$n==0}
 denoms<-z[!duplicated(z$site),]
 den<-if(safe)sum(denoms$denominator_n) else NA_real_;obs<-if(safe)sum(denoms$observed_n) else NA_real_;miss<-if(safe)sum(denoms$missing_n) else NA_real_
 for(category in unique(z$category)){
  q<-z[z$category==category,];count<-if(safe)sum(q$n) else NA_real_
  categorical[[paste(k,category)]]<-tibble::tibble(cohort=first$cohort,group=first$group,variable=first$variable,category=category,n=count,denominator_n=den,observed_n=obs,missing_n=miss,percent_total=if(safe&&den>0)100*count/den else NA_real_,status=if(safe)'observed_all_sites' else'incomplete_or_suppressed_site_partition')
 }
}
pool_cat<-dplyr::bind_rows(categorical);readr::write_csv(pool_cat,file.path(out,'pooled_categorical.csv'),na='')
quant<-list()
for(k in unique(do.call(paste,c(all$histograms[keys],sep='|')))){
 z<-all$histograms[do.call(paste,c(all$histograms[keys],sep='|'))==k,];first<-z[1,]
 safe<-expected(first$cohort,first$group,first$variable)&&length(unique(z$site))==length(sites)&&all(z$cell_status=='observed')
 if(safe){h<-z |> dplyr::group_by(bin,lower,upper,integer_exact) |> dplyr::summarise(n=sum(n),.groups='drop') |> dplyr::arrange(bin);N<-sum(h$n)}else N<-0
 for(p in c(.25,.5,.75)){
  estimate<-lo<-hi<-NA_real_;status<-'incomplete_or_suppressed_histogram'
  if(safe&&N>0){r<-1+(N-1)*p;r1<-floor(r);r2<-ceiling(r);j1<-which(cumsum(h$n)>=r1)[1];j2<-which(cumsum(h$n)>=r2)[1];lo<-h$lower[j1];hi<-h$upper[j2]
   if(all(h$integer_exact)&&all(is.finite(c(lo,hi)))){v1<-(h$lower[j1]+h$upper[j1])/2;v2<-(h$lower[j2]+h$upper[j2])/2;estimate<-v1+(r-r1)*(v2-v1);status<-'exact_from_integer_counts'}else status<-'fixed_bin_bounds_only'
  }else if(safe)status<-'all_missing'
  quant[[paste(k,p)]]<-tibble::tibble(cohort=first$cohort,group=first$group,variable=first$variable,probability=p,estimate=estimate,lower_bound=lo,upper_bound=hi,status=status)
 }
}
readr::write_csv(dplyr::bind_rows(quant),file.path(out,'pooled_quantiles.csv'),na='')
for(k in c('continuous','categorical')){
 z<-all[[k]];readr::write_csv(z[z$cohort%in%c('eligible_cohort','landmark_cohort')&z$group=='Overall',],file.path(out,paste0('table1_',k,'_by_site.csv')),na='')
 readr::write_csv(z[z$cohort=='primary_phenotypes',],file.path(out,paste0('table2_',k,'_by_site.csv')),na='')
}
# Reconcile every released pooled categorical partition.
if(nrow(pool_cat)){q<-pool_cat |> dplyr::filter(status=='observed_all_sites') |> dplyr::group_by(cohort,group,variable) |> dplyr::summarise(total=sum(n),den=dplyr::first(denominator_n),.groups='drop');stopifnot(all(q$total==q$den))}
jsonlite::write_json(list(status='passed',sites=sites,consensus_id=m[[1]]$consensus_id,patient_level_inputs_read=FALSE,mean_SD_pooling='exact observed n,sum,sum_squares',categorical_denominators_counted_once_per_site=TRUE,site_medians_not_averaged=TRUE,suppression_not_imputed=TRUE,patient_counts_not_assumed_globally_unique=TRUE,quantiles='integer-count exact where released; continuous fixed-bin bounds only; suppressed histograms unavailable'),file.path(out,'collection_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(c('Publication descriptor return: Table1/Table2','',
'Table1: eligible-cohort and primary-landmark overall characteristics by site; keep these cohorts separate.',
'Table2: frozen shared primary P groups plus I and Overall; do not use provisional local-cluster tables as shared phenotypes.',
'Files include demographics/year, noncancer Charlson and coded domains, available oncology proxies/fields, original modified SOFA, adapted APACHE proxy with missingness, support hours/any-use/last-hour state, ARF evidence and unknown hours, code status/transitions, whole-stay descriptive outcomes, and separate28d post-landmark outcomes. Hospital LOS is descriptive only.',
'Exact observed n/sum/sum-squares support pooled means/SD; categorical n/explicit denominators support pooled percentages. Missing/unavailable/suppressed status is explicit. Never average site medians or report a suppressed missing cell as zero.',
'Pooled medians/IQR require released common-bin counts: exact for integer-valued variables when possible; otherwise report bounds or use site-specific medians/IQR. Do not turn coarse bins into exact medians.',
'Unique patient counts are site-distinct and not additive across phenotype groups or institutions; hospitalization counts are the analysis denominator.',
'Other site confirmations still needed: IRB/waiver/disclosures, ETL/timezone/provenance, treatment-limit mapping, and optional oncology availability/timing. Unavailable stage/ECOG/therapy data and validated APACHE IV cannot be filled by these descriptors.'),file.path(out,'descriptor_return_report.txt'))
cat('Table1/Table2 aggregate collection and reconciliation passed\n')
