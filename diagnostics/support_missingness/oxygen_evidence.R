source("code/lib/load.R")
out <- Sys.getenv("CLIF_SUPPORT_AUDIT_OUTPUT", unset = "output/run_UCMC_support_audit_20261022")
h <- readRDS(file.path(out,"hour_audit_local.rds"))
d <- readRDS(file.path(init_run(read_run_config()),"local","derivation.rds"))$analyses$h72_primary$data
r <- readRDS(file.path(out,"respiratory_source_local.rds"))
r$state <- support_state(r$device_category)
rb <- split(r,r$hospitalization_id)
checks <- lapply(seq_len(nrow(d)),function(i) {
 id <- d$hospitalization_id[i]; z <- h[h$hospitalization_id==id & h$support=="UNKNOWN",]
 if(!nrow(z))return(NULL)
 x <- rb[[id]]; if(is.null(x))x<-r[0,]
 x <- x[which(x$recorded_dttm>=d$admission_dttm[i] & x$recorded_dttm<d$landmark[i]),]
 x <- x[order(x$recorded_dttm),]
 dplyr::bind_rows(lapply(seq_len(nrow(z)),function(j) {
  start <- d$t0[i]+z$h[j]*3600; end <- start+3600
  now <- x[which(x$recorded_dttm>=start & x$recorded_dttm<end),]
  valid_now <- now[which(!is.na(now$lpm_set) & now$lpm_set>=0 & now$lpm_set<=60),]
  observed <- any(valid_now$lpm_set>0)
  past <- x[which(x$recorded_dttm<=start & x$recorded_dttm>=start-6*3600 &
   !is.na(x$lpm_set) & x$lpm_set>=0 & x$lpm_set<=60),]
  last_named <- x[which(x$recorded_dttm<=start & x$state!="UNKNOWN"),]
  carried <- FALSE
  if(!nrow(valid_now) && nrow(past)) {
   last <- tail(past,1)
   carried <- last$lpm_set>0 && (!nrow(last_named) || last$recorded_dttm>=tail(last_named$recorded_dttm,1))
  }
  tibble::tibble(hospitalization_id=id,h=z$h[j],local_cluster=z$local_cluster[j],
   positive_flow_observed=observed,positive_flow_carried=carried)
 }))
})
checks <- dplyr::bind_rows(checks)
saveRDS(checks,file.path(out,"oxygen_evidence_local.rds"))
summary <- checks |> dplyr::group_by(local_cluster) |>
 dplyr::summarise(unknown_hours=dplyr::n(),observed_oxygen_hours=sum(positive_flow_observed),
  carried_oxygen_hours=sum(!positive_flow_observed & positive_flow_carried),
  residual_unknown_hours=sum(!positive_flow_observed & !positive_flow_carried),.groups="drop")
print(summary,width=Inf);readr::write_csv(summary,file.path(out,"oxygen_evidence_summary.csv"))
cat("Upper bound diagnostic only: oxygen evidence does not identify support modality. Primary states unchanged.\n")
stopifnot(nrow(h)==nrow(d)*72,sum(h$support=="UNKNOWN")==sum(d$hours_UNKNOWN),nrow(checks)==sum(d$hours_UNKNOWN))
