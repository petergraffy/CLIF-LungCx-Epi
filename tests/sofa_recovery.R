source('code/lib/load.R');source('diagnostics/sofa_completeness/recovery.R')
time<-utc('2020-01-01 00:00:00')
v<-tibble::tibble(hospitalization_id=c('a','a','b','b','b','c','c','d','d'),recorded_dttm=time+c(0,0,0,0,0,0,60,0,0),vital_category=c('sbp','dbp','sbp','sbp','dbp','sbp','dbp','sbp','dbp'),vital_value=c(120,60,120,130,60,120,60,60,100))
g<-tibble::tibble(hospitalization_id=rep(c('a','b','c'),each=3),recorded_dttm=time,assessment_category=rep(c('gcs_eye','gcs_verbal','gcs_motor'),3),numerical_value=c(4,5,6,4,NA,6,4,5,7))
x<-sofa_recover_tuples(list(vitals=v,patient_assessments=g))
stopifnot(nrow(x$map)==1,x$map$hospitalization_id=='a',x$map$vital_value==80,nrow(x$gcs)==1,x$gcs$hospitalization_id=='a',x$gcs$numerical_value==15)
# Conflicting and asynchronous components must not yield a calculated score.
cat('SOFA recovery measurement contracts passed\n')
