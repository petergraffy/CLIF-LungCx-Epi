# Aggregate Aalen-Johansen exports on a fixed grid; no event timestamps returned.
released_cif_grid<-function(d,days=0:28,min_cell=10L){
 min_cell<-aggregate_min_cell(min_cell)
 stopifnot(!anyNA(d$phenotype),all(d$time_from_landmark>0))
 cif<-aj_grid(d,days=days)
 cif$cell_status<-'observed';cif$estimate<-cif$cumulative_incidence;cif$cumulative_incidence<-NULL
 for(ph in unique(cif$phenotype))for(day in unique(cif$day[cif$phenotype==ph])){
  ix<-which(cif$phenotype==ph & cif$day==day);z<-d[d$phenotype==ph,]
  counts<-table(factor(z$event_type[z$time_from_landmark<=day],levels=c('death','hospice','live_discharge')))
  unsafe<-any(counts>0 & counts<min_cell) || any(cif$n_at_risk[ix]>0 & cif$n_at_risk[ix]<min_cell)
  if(unsafe){cif$estimate[ix]<-NA_real_;cif$std_error[ix]<-NA_real_;cif$cell_status[ix]<-'suppressed_partition'}
 }
 cif$n<-cif$denominator;cif$n_at_risk_exact<-cif$n_at_risk
 band<-t(vapply(cif$denominator,aggregate_count_band,numeric(2)));cif$n_lower<-band[,1];cif$n_upper<-band[,2]
 band<-t(vapply(cif$n_at_risk,aggregate_count_band,numeric(2)));cif$at_risk_lower<-band[,1];cif$at_risk_upper<-band[,2]
 cif$denominator<-NULL;cif$n_at_risk<-NULL
 cif
}
