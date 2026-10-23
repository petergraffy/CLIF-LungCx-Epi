source("diagnostics/federated/model_lib.R")
site<-commandArgs(trailingOnly=TRUE)[1];root<-fed_root();home<-file.path(root,"sites",site)
d<-readRDS(file.path(home,"local","full","outcomes.rds"))$data;rows<-list()
for(ph in sort(unique(d$phenotype))){z<-d[d$phenotype==ph,];n<-nrow(z);nb<-fed_band_n(n)
 yd<-as.integer(z$time_from_landmark<=28 & z$event_type=="death");yh<-as.integer(z$time_from_landmark<=28 & z$event_type=="hospice")
 parts<-c(sum(yd),sum(yh),n-sum(yd)-sum(yh));cause_safe<-!any(parts>0 & parts<10)
 for(endpoint in c("death","hospice","death_or_hospice")){
 y<-switch(endpoint,death=yd,hospice=yh,death_or_hospice=yd+yh);r<-mean(y)
 safe<-n>=20 && sum(y)>=10 && sum(1-y)>=10 && (endpoint=="death_or_hospice" || cause_safe)
 score<-tapply(y-r,z$patient_id,sum);g<-length(score);se<-if(g>1)sqrt(sum(score^2)*g/(g-1))/n else NA_real_
 rows[[paste(ph,endpoint)]]<-tibble::tibble(site=site,phenotype=ph,endpoint=endpoint,n_lower=nb[1],n_upper=nb[2],
 risk=if(safe)r else NA_real_,std_error=if(safe)se else NA_real_,conf_low=if(safe)max(0,r-1.96*se) else NA_real_,conf_high=if(safe)min(1,r+1.96*se) else NA_real_,status=if(safe)"ok" else "withheld_small_partition")
 }}
fed_write(dplyr::bind_rows(rows),file.path(home,"export","risks","risks_28day.csv"))
message(site," local 28-day absolute risks exported with rare partitions suppressed")
