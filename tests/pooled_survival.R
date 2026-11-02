source('diagnostics/pooled_survival/lib.R')
T1<-c(1,1,2);E1<-c(1L,0L,1L);T2<-c(1.5,3,4,5);E2<-c(0L,1L,0L,1L)
grid<-c(0,1,2,3,5,7)
a<-ps_counts(T1,E1,grid);b<-ps_counts(T2,E2,grid)
pooled<-ps_product(a$risk+b$risk,a$events+b$events)
fit<-survival::survfit(survival::Surv(c(T1,T2),c(E1,E2))~1,timefix=FALSE)
expected<-summary(fit,times=grid,extend=TRUE)$surv
stopifnot(all(abs(pooled-expected)<1e-12))
# Unequal sites and different censoring: averaging site KM is not pooled KM.
average<-3/7*ps_product(a$risk,a$events)+4/7*ps_product(b$risk,b$events)
stopifnot(abs(average[4]-pooled[4])>.01)
# Bootstrap weights keep both admissions for the same patient together.
w1<-c(2L,2L,1L);w2<-c(1L,2L,1L,1L)
a<-ps_counts(T1,E1,grid,w1);b<-ps_counts(T2,E2,grid,w2)
fit<-survival::survfit(survival::Surv(c(T1,T2),c(E1,E2))~1,weights=c(w1,w2),timefix=FALSE)
stopifnot(all(abs(ps_product(a$risk+b$risk,a$events+b$events)-summary(fit,times=grid,extend=TRUE)$surv)<1e-12))
empty<-ps_counts(numeric(),integer(),grid);stopifnot(all(empty$risk==0),all(empty$events==0),all(is.na(ps_product(empty$risk,empty$events))))
# CSV decimal round trips must not lose tied events or shift their risk sets.
T<-ps_time(c(1234567,1234567,2345678));g<-c(0,T[c(1,3)],28)
roundtrip<-as.numeric(sprintf('%.16g',g));q<-ps_counts(T,c(1L,0L,1L),roundtrip)
stopifnot(identical(q$risk,c(3,3,1,0)),identical(q$events,c(0,1,1,0)))
cat('PASS: federated risk/event pooling equals centralized KM; event/censor ties, unequal sites, weighted patient resamples and empty sites\n')
# A death prevents later live discharge; censoring it would overestimate incidence.
d<-data.frame(time=c(1,1,2,3,4),event_type=c('death','live_discharge','hospice','live_discharge','censored'),valid_outcome=c(TRUE,TRUE,TRUE,TRUE,FALSE))
g<-c(0,1,2,3,4,7)
live<-ps_counts(d$time,ps_event(d,'live_discharge'),g)
all_exit<-ps_counts(d$time,ps_event(d,'any_exit'),g)
ci<-ps_incidence(live$risk,live$events,all_exit$events)
fit<-survival::survfit(survival::Surv(time,factor(event_type,levels=c('censored','death','hospice','live_discharge')))~1,data=d,timefix=FALSE)
expected<-summary(fit,times=g,extend=TRUE)$pstate[,'live_discharge']
stopifnot(all(abs(ci-expected)<1e-12),abs(tail(ci,1)-.4)<1e-12)
stopifnot(tail(1-ps_product(live$risk,live$events),1)>tail(ci,1))
stopifnot(all(is.na(ps_incidence(rep(0,3),rep(0,3),rep(0,3)))))
cat('PASS: discharge incidence matches multistate Aalen-Johansen with tied death/live exits, hospice and censoring\n')

# Composite death/hospice matches multistate AJ and complements live discharge.
adverse<-ps_counts(d$time,ps_event(d,'death_or_hospice'),g)
adverse_ci<-ps_incidence(adverse$risk,adverse$events,all_exit$events)
expected_adverse<-rowSums(summary(fit,times=g,extend=TRUE)$pstate[,c('death','hospice')])
stopifnot(all(abs(adverse_ci-expected_adverse)<1e-12))
stopifnot(all(abs(adverse_ci+ci+ps_product(all_exit$risk,all_exit$events)-1)<1e-12))
cat('PASS: composite death/hospice incidence matches multistate AJ; both causes plus event-free probability sum to one\n')
