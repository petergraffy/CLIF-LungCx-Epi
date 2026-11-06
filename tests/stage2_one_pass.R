source('scripts/deployment_lib.R');source('diagnostics/stage2/lib.R');source('diagnostics/pooled_survival/lib.R')
ids<-rep(paste0('synthetic_',1:40),each=2);fold<-s2_folds(ids)
stopifnot(all(vapply(split(fold,ids),function(x)length(unique(x))==1L,logical(1))))
x<-cbind(`(Intercept)`=1,signal=rep(seq(-1,1,length.out=40),each=2));y<-rep(rep(c(0,1),20),each=2)
f<-s2_oof(x,y,ids,fold);stopifnot(all(is.finite(f$prediction)),all(f$audit$patient_overlap==0),all(f$audit$status=='ok'))
# No-outcome training must never produce apparent or fabricated validation scores.
f<-s2_oof(x,rep(0,length(y)),ids,fold);stopifnot(all(is.na(f$prediction)),all(f$audit$status=='nonestimable_training_fold'))
# Histograms reconstruct true union-time risk sets, including tied censors/events.
times<-c(1,1,2,3,4,28);grid<-c(0,1,1.5,2,3,4,7,28);weights<-c(1,2,1,0,3,2)
r<-s2_hist_risk(times,weights,grid,sum(weights));expected<-vapply(grid,function(t)sum(weights[times>=t]),numeric(1));stopifnot(identical(as.numeric(r),expected))
event<-c(1,0,1,1,0,0);q<-ps_counts(times,event,grid,weights);stopifnot(identical(r,q$risk))
# CSV templates must reconstruct exactly with arbitrary serialized feature order.
b<-tempfile(tmpdir='output');dir.create(b);p<-matrix(c(.2,.8,.9,.1),2,dimnames=list(c('P1','P2'),c('a','b')))
z<-expand.grid(phenotype=rownames(p),feature=colnames(p));z$probability<-p[cbind(match(z$phenotype,rownames(p)),match(z$feature,colnames(p)))];readr::write_csv(z,file.path(b,'prototypes.csv'));stopifnot(isTRUE(all.equal(deploy_load_prototypes(b),p)))
cat('PASS: one-pass grouped OOF leakage checks, nonestimable folds, exact tied-time aggregate risk reconstruction, portable CSV prototypes\n')
