source('diagnostics/revision20261022/lib.R')
site<-commandArgs(trailingOnly=TRUE)[1];old<-'output/run_federated_rehearsal_20261022';new<-old_root();root<-'output/run_three_site_refresh_20261026'
a<-readRDS(file.path(old,'sites',site,'local/discovery.rds'));b<-readRDS(file.path(new,'sites',site,'local/discovery.rds'))
x<-readRDS(file.path(old,'sites',site,'local/full/outcomes.rds'))$data;y<-readRDS(file.path(new,'sites',site,'local/full/outcomes.rds'))$data
stopifnot(setequal(x$hospitalization_id,y$hospitalization_id));y<-y[match(x$hospitalization_id,y$hospitalization_id),]
if(site!='UCMC')stopifnot(identical(a$panel,b$panel))
rows<-lapply(c('previous','refreshed'),function(version){d<-if(version=='previous')x else y;p<-if(version=='previous')a$panel else b$panel
 tibble::tibble(site=site,version=version,n=nrow(d),patients=dplyr::n_distinct(d$patient_id),unknown_hours=sum(p$support=='UNKNOWN'),total_hours=nrow(p),unknown_pct=100*mean(p$support=='UNKNOWN'),arf_positive_n=sum(d$measured_arf),complete_sofa_n=sum(!is.na(d$sofa_total)),local_k=if(version=='previous')a$fit$k else b$fit$k,classified_n_lower=fed_band_n(sum(d$phenotype!='I'))[1],classified_n_upper=fed_band_n(sum(d$phenotype!='I'))[2])})
out<-file.path(root,'sites',site,'export');fed_write(dplyr::bind_rows(rows),file.path(out,'source_comparison.csv'))
cross<-tibble::tibble(previous=x$phenotype,refreshed=y$phenotype) |> dplyr::count(previous,refreshed,name='n') |> dplyr::mutate(n_lower=if(aggregate_exact())n else 10*floor(n/10),n_upper=if(aggregate_exact())n else 10*floor(n/10)+9) |> dplyr::select(-n)
fed_write(cross |> dplyr::mutate(site=site),file.path(out,'assignment_transitions.csv'))
ari<-adjusted_rand(as.integer(factor(x$phenotype)),as.integer(factor(y$phenotype)));ix<-x$phenotype!='I'&y$phenotype!='I'
ari2<-adjusted_rand(as.integer(factor(x$phenotype[ix])),as.integer(factor(y$phenotype[ix])))
fed_write(tibble::tibble(site=site,adjusted_rand_all=ari,adjusted_rand_both_classified=ari2,interpretation='Agreement of independently rederived partitions, not clinical accuracy or necessarily equivalent labels'),file.path(out,'partition_comparison.csv'))
fed_json(list(site=site,status='passed',same_landmark_members=TRUE,source_panel_unchanged=identical(a$panel,b$panel),all_refreshed_cases_retained=TRUE,new_input=b$input_root),file.path(out,'comparison_audit.json'))
