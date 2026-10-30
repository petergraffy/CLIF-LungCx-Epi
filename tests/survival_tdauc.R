source('code/lib/load.R');source('diagnostics/manuscript_figures/survival_tdauc_lib.R')
T<-c(rep(2,20),rep(4,20),rep(40,20));delta<-c(rep(1L,20),rep(2L,20),rep(1L,20));p<-c(rep(.9,20),rep(.5,20),rep(.1,20))
stopifnot(td_auc_rank(T,delta,p,7,1)==1,td_auc_rank(T,delta,p,7,2)==1)
stopifnot(td_auc_rank(T,delta,rep(.5,60),7,1)==.5,td_auc_rank(T,delta,rep(.5,60),7,2)==.5)
# Competing discharges count as controls only in definition2.
p2<-c(rep(.6,20),rep(.9,20),rep(.1,20))
stopifnot(td_auc_rank(T,delta,p2,7,1)==1,td_auc_rank(T,delta,p2,7,2)==.5)
q<-timeROC::timeROC(T=T,delta=delta,marker=p2,cause=1,times=c(7,14,28)+1e-8,iid=FALSE)
stopifnot(all(abs(q$AUC_1-1)<1e-12),all(abs(q$AUC_2-.5)<1e-12),all(q$weights$IPCW.times==1))
# All hospital exits observed before a horizon: AUC_2 remains identified, whereas
# event-free-control AUC_1 has no controls and must remain unavailable.
q<-timeROC::timeROC(T=c(rep(2,20),rep(4,20)),delta=c(rep(1,20),rep(2,20)),marker=c(rep(.9,20),rep(.1,20)),cause=1,times=c(3,7),iid=FALSE)
stopifnot(is.na(q$AUC_1[2]),q$AUC_2[2]==1)
d<-tibble::tibble(hospitalization_id=paste0('h',1:60),patient_id=paste0('p',1:60),valid_outcome=TRUE,event_type=ifelse(delta==1,'death','live_discharge'))
z<-km_release(d,T)
stopifnot(all(z$survival[z$day==0]==1),all(z$at_risk_lower[z$day==0]==60),all(z$survival[z$day==7]==2/3))
d$valid_outcome[1]<-FALSE
z<-km_release(d,T)
stopifnot(abs(z$survival[z$day==7&z$endpoint=='death']-41/60)<1e-12)
cat('PASS: KM starts at100%, censoring removes risk without a death, tied-marker tdAUC, competing-control definitions and official timeROC agreement\n')
