source('diagnostics/revision20261022/lib.R')
site<-commandArgs(trailingOnly=TRUE)[1];home<-file.path(rev_root(),'sites',site);d<-readRDS(file.path(home,'local/review.rds'))$data
vars<-c('age_years','charlson_noncancer','sofa_total','hours_IMV','hours_HFNC','hours_CPAP','hours_NIPPV','hours_UNKNOWN','hypox_fraction','hyper_fraction','mixed_fraction','positive_fraction');rows<-list()
for(ph in unique(d$phenotype))for(v in vars){
 z<-d[d$phenotype==ph,v,drop=TRUE];nn<-sum(!is.na(z));nb<-fed_band_n(nn);mb<-fed_band_n(sum(is.na(z)));safe<-nn>=max(1,aggregate_min_cell(10))
 rows[[paste(ph,v)]]<-tibble::tibble(site=site,phenotype=ph,variable=v,observed_lower=nb[1],observed_upper=nb[2],missing_lower=mb[1],missing_upper=mb[2],mean=if(safe)mean(z,na.rm=TRUE) else NA_real_,median=if(safe)median(z,na.rm=TRUE) else NA_real_,q1=if(safe)unname(quantile(z,.25,na.rm=TRUE)) else NA_real_,q3=if(safe)unname(quantile(z,.75,na.rm=TRUE)) else NA_real_,status=if(safe)'ok' else 'withheld_small_observed_partition')
}
fed_write(dplyr::bind_rows(rows),file.path(home,'export/continuous_denominators.csv'))
d$outcome28<-ifelse(d$time_from_landmark<=28,d$event_type,'no_terminal_event_by_day28')
fed_write(band_table(d,c('phenotype','code_at_landmark'),'outcome28',d$outcome28) |> dplyr::mutate(site=site),file.path(home,'export/outcomes_by_code.csv'))

# Apply the observed-denominator release rule to the wide convenience table too.
f<-file.path(home,'export/clinical_table.csv');ct<-readr::read_csv(f,show_col_types=FALSE)
for(i in seq_len(nrow(ct)))for(v in intersect(names(ct),vars))if(sum(!is.na(d[[v]][d$phenotype==ct$phenotype[i]]))<aggregate_min_cell(10))ct[[v]][i]<-NA_real_
fed_write(ct,f)
