options(clif.aggregate_policy='legacy') # Historical export compatibility test
source('code/lib/load.R');source('diagnostics/manuscript_figures/cif_lib.R')
# With complete follow-up and no censoring, AJ equals event proportions.
n<-60L
d<-tibble::tibble(hospitalization_id=paste0('h',seq_len(n)),patient_id=paste0('p',seq_len(n)),phenotype='P1',time_from_landmark=rep(c(2,4,6),each=20),event_type=rep(c('death','hospice','live_discharge'),each=20))
x<-released_cif_grid(d,0:28,10)
stopifnot(max(x$day)==6,all(x$cell_status=='observed'),!any(c('patient_id','hospitalization_id','time_from_landmark')%in%names(x)))
z<-x[x$day==6,];stopifnot(abs(sum(z$estimate)-1)<1e-12,all(abs(z$estimate[z$event!='(s0)']-1/3)<1e-12))
# A rare competing event must mask all states at that day, not only its own CIF.
d$event_type[1]<-'hospice';x<-released_cif_grid(d,0:6,10)
stopifnot(all(x$cell_status[x$day==2]=='suppressed_partition'),all(is.na(x$estimate[x$day==2])),all(is.na(x$std_error[x$day==2])))
stopifnot(all(x$cell_status[x$day==4]=='observed'))
cat('PASS: competing-event probabilities, linked suppression, fixed grid, no extrapolation and no identifiers\n')
