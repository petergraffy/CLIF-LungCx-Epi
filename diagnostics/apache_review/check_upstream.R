# Read/review a pinned helper; never invoke its cohort-calculation function.
# Source snapshots stay in ignored output, not redistributed in application code.
suppressPackageStartupMessages(library(data.table))
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2)stop('Usage: check_upstream.R REVIEWED_HELPER NEW_REPORT_DIRECTORY')
path<-args[1];out<-args[2];dir.create(out,recursive=TRUE,showWarnings=FALSE)
x<-parse(path)
stopifnot(all(vapply(x,function(z)is.call(z)&&as.character(z[[1]])%in%c('=','<-')&&is.symbol(z[[2]])&&is.call(z[[3]])&&as.character(z[[3]][[1]])=='function',logical(1))))
# This reviewed helper contains only point-function definitions and calculate_apache2.
expected<-c('pts_temp','pts_map','pts_hr','pts_rr','pts_ph','pts_sodium','pts_potassium','pts_creatinine','pts_hemoglobin','pts_wbc','pts_oxygenation','pts_gcs','pts_age','calculate_apache2')
stopifnot(setequal(vapply(x,function(z)as.character(z[[2]]),character(1)),expected))
e<-new.env(parent=asNamespace('data.table'));for(z in x)eval(z,e)
# Static definition loading only; do not run setup, load clinical data or call
# calculate_apache2. These probes document behavior, not clinical validation.
z<-data.frame(check=c('hemoglobin_10_as_supplied_g_dL','hemoglobin_100_after_times10','temperature_38_45','MAP_109_5','sodium_129_5','GCS_missing','creatinine4_default','creatinine4_explicit_ARF'),points=c(e$pts_hemoglobin(10),e$pts_hemoglobin(100),e$pts_temp(38.45),e$pts_map(109.5),e$pts_sodium(129.5),e$pts_gcs(NA),e$pts_creatinine(4),e$pts_creatinine(4,TRUE)))
write.csv(z,file.path(out,'point_function_probes.csv'),row.names=FALSE,na='')
resp<-data.table(id='a',source_time=as.POSIXct('2020-01-01 10:00:00',tz='UTC'),fio=.8)
abg<-data.table(id='a',test_time=as.POSIXct(c('2020-01-01 09:00:00','2020-01-01 00:00:00','2020-01-01 15:00:00'),tz='UTC'))
setkey(resp,id,source_time)
setnames(abg,'test_time','source_time')
joined<-resp[abg,roll=4*3600,rollends=c(TRUE,TRUE)]
write.csv(data.frame(check=c('ABG_one_hour_before_first_FiO2','ABG_ten_hours_before','ABG_five_hours_after'),matched_fio2=joined$fio),file.path(out,'rolling_join_probes.csv'),row.names=FALSE,na='')
stopifnot(z$points[1]==4,z$points[2]==0,all(is.na(z$points[3:6])),z$points[7]==4,z$points[8]==8,joined$fio[1]==.8,is.na(joined$fio[2]))
writeLines(c('Reviewed helper behavior reproduced. These findings require correction or explicit proxy assumptions before integration.','The full cohort function was not run and no patient records were loaded.','Data.table version:',as.character(packageVersion('data.table')),'Decoded source snapshot MD5:',unname(tools::md5sum(path))),file.path(out,'probe_status.txt'))
print(z)
