# One-exchange Stage 2 primitives. All patient rows stay at the site.
source('diagnostics/apache_adapted/lib.R')
s2_folds<-function(patient_id,k=5L,seed=20261105L){
 ids<-sort(unique(patient_id),method='radix');k<-min(k,length(ids))
 if(k<2L)return(rep(NA_integer_,length(patient_id)))
 set.seed(seed);fold<-sample(rep(seq_len(k),length.out=length(ids)))
 fold[match(patient_id,ids)]
}
s2_oof<-function(x,y,patient_id,fold){
 p<-rep(NA_real_,length(y));audit<-list()
 for(f in sort(unique(fold[!is.na(fold)]))){
  test<-which(fold==f);train<-which(fold!=f)
  stopifnot(!any(patient_id[test]%in%patient_id[train]))
  fit<-if(length(unique(y[train]))==2L)tryCatch(ridge_logistic(x[train,,drop=FALSE],y[train],patient_id[train]),error=identity) else NULL
  good<-is.list(fit)&&!inherits(fit,'error')&&isTRUE(fit$converged)
  if(good)p[test]<-plogis(drop(x[test,,drop=FALSE]%*%fit$beta))
  audit[[as.character(f)]]<-tibble::tibble(fold=f,training_n=length(train),validation_n=length(test),training_patients=length(unique(patient_id[train])),validation_patients=length(unique(patient_id[test])),training_events=sum(y[train]),status=if(good)'ok' else 'nonestimable_training_fold',patient_overlap=0L)
 }
 list(prediction=p,audit=dplyr::bind_rows(audit))
}
# Aggregate exit histogram -> exact risk/event counts at any union of site times.
# Censors tied with an event remain in its risk set.
s2_hist_risk<-function(time,count,grid,n){
 ord<-order(time);idx<-findInterval(grid,time[ord],left.open=TRUE)
 n-c(0,cumsum(count[ord]))[idx+1L]
}
