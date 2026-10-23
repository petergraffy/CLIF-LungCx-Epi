source("code/lib/load.R")
root <- init_run(read_run_config())
out <- Sys.getenv("CLIF_SUPPORT_AUDIT_OUTPUT", unset = "output/run_UCMC_support_audit_20261022")
p <- readRDS(file.path(root,"local","prepared.rds"))
a <- readRDS(file.path(root,"local","derivation.rds"))$analyses$h72_primary
rs <- readRDS(file.path(out,"respiratory_source_local.rds"))
rs$state <- support_state(rs$device_category)
rb <- split(rs,rs$hospitalization_id)
ab <- split(p$tables$adt,p$tables$adt$hospitalization_id)
vb <- split(p$tables$vitals,p$tables$vitals$hospitalization_id)
rows <- lapply(seq_len(nrow(a$data)),function(i) {
 d <- a$data[i,]; id <- d$hospitalization_id
 r <- rb[[id]]; if(is.null(r))r<-rs[0,]; ad <- ab[[id]]
 known <- r[r$state!="UNKNOWN" & !is.na(r$recorded_dttm),]
 known <- known[order(known$recorded_dttm),]
 sq <- a$fit$sequence[match(id,a$fit$ids),]
 support <- vapply(strsplit(sq,"|",fixed=TRUE),function(x)x[1],character(1))
 dplyr::bind_rows(lapply(0:71,function(h) {
  start <- d$t0+h*3600; end <- start+3600
  inhour <- r[which(r$recorded_dttm>=start & r$recorded_dttm<end),]
  before <- known[which(known$recorded_dttm<=start),]
  after <- known[which(known$recorded_dttm>=end),]
  prev_age <- if(nrow(before)) hours_between(start,tail(before$recorded_dttm,1)) else Inf
  next_age <- if(nrow(after)) hours_between(after$recorded_dttm[1],end) else Inf
  loc <- ad[which(ad$in_dttm<end & (is.na(ad$out_dttm)|ad$out_dttm>start)),]
  location <- if(any(tolower(loc$location_category)=="icu"))"icu" else
   if(any(tolower(loc$location_category)=="ward"))"ward" else "other_or_no_adt"
  cause <- if(support[h+1]!="UNKNOWN") "known_support" else if(nrow(inhour))
   "source_rows_without_device" else "no_respiratory_row_in_hour"
  gap <- if(support[h+1]!="UNKNOWN") "known_support" else if(!is.finite(prev_age))
   "before_first_known_device" else if(prev_age>6)"more_than_6h_since_known_device" else "other"
  tibble::tibble(hospitalization_id=id,local_cluster=d$local_cluster,h=h,
   location=location,support=support[h+1],cause=cause,gap=gap,
   previous_known_hours=prev_age,next_known_hours=next_age,
   fio2_present=any(!is.na(inhour$fio2_set)),positive_flow=any(inhour$lpm_set>0,na.rm=TRUE),
   zero_flow=any(inhour$lpm_set==0,na.rm=TRUE),peep_present=any(!is.na(inhour$peep_set)),
   airway_present=any(!is.na(inhour$artificial_airway)),
   bracketing_same_device=nrow(before)>0 && nrow(after)>0 && tail(before$state,1)==after$state[1])
 }))
})
hours <- dplyr::bind_rows(rows)
saveRDS(hours,file.path(out,"hour_audit_local.rds"))
for (group in list(c("local_cluster","support"),c("location","support"),
 c("cause","location"),c("gap","location"))) {
 x <- hours
 if(any(c("cause","gap")%in%group))x<-x[x$support=="UNKNOWN",]
 z<-x |> dplyr::count(dplyr::across(dplyr::all_of(group)),name="hours")
 readr::write_csv(z,file.path(out,paste0(paste(group,collapse="_"),".csv")))
 print(z,n=50)
}
z<-hours |> dplyr::filter(support=="UNKNOWN") |> dplyr::group_by(location) |>
 dplyr::summarise(hours=dplyr::n(),fio2_hours=sum(fio2_present),positive_flow_hours=sum(positive_flow),
 zero_flow_hours=sum(zero_flow),peep_hours=sum(peep_present),airway_hours=sum(airway_present),
 same_device_bracketed_hours=sum(bracketing_same_device),.groups="drop")
print(z,width=Inf);readr::write_csv(z,file.path(out,"unknown_source_evidence.csv"))
print(hours |> dplyr::group_by(hospitalization_id,local_cluster) |>
 dplyr::summarise(unknown_hours=sum(support=="UNKNOWN"),.groups="drop") |>
 dplyr::group_by(local_cluster) |> dplyr::summarise(n=dplyr::n(),entirely_unknown=sum(unknown_hours==72),
 median_unknown=median(unknown_hours),.groups="drop"))
cat("Vital categories:\n");print(p$tables$vitals |> dplyr::count(vital_category,sort=TRUE),n=50)
