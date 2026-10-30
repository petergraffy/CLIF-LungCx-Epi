# Public score bins; no individual scores, identifiers or event times exported.
ROC_BINS<-100000L
roc_bin<-function(p){stopifnot(all(is.finite(p)),all(p>=0&p<=1));pmin(ROC_BINS-1L,floor(p*ROC_BINS))}
roc_hist<-function(bin,event,weight=rep(1,length(bin)),grid=sort(unique(bin))){
 stopifnot(length(bin)==length(event),length(weight)==length(bin),!anyNA(match(bin,grid)))
 sum_to_grid<-function(w){out<-numeric(length(grid));if(length(w)){q<-rowsum(w,match(bin,grid),reorder=FALSE);out[as.integer(rownames(q))]<-q[,1]};out}
 list(bin=grid,cases=sum_to_grid(weight*event),controls=sum_to_grid(weight*(1-event)))
}
roc_auc<-function(cases,controls){
 if(!sum(cases)||!sum(controls))return(NA_real_)
 sum(cases*(cumsum(controls)-controls+.5*controls))/(sum(cases)*sum(controls))
}
