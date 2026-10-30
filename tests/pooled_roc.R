source('diagnostics/pooled_roc/lib.R')
p<-c(.1,.1,.9,.6,.2,.6,.8,1);y<-c(1,0,1,0,0,1,0,1)
grid<-sort(unique(roc_bin(p)));a<-roc_hist(roc_bin(p[1:3]),y[1:3],grid=grid);b<-roc_hist(roc_bin(p[4:8]),y[4:8],grid=grid)
auc<-roc_auc(a$cases+b$cases,a$controls+b$controls)
expected<-(sum(rank(p)[y==1])-sum(y)*(sum(y)+1)/2)/(sum(y)*sum(1-y))
stopifnot(abs(auc-expected)<1e-12,roc_bin(1)==99999,roc_bin(0)==0)
w<-c(2,2,1,1,1,1,2,1);h<-roc_hist(roc_bin(p),y,w)
idx<-rep(seq_along(p),w);z<-p[idx];v<-y[idx]
e<-(sum(rank(z)[v==1])-sum(v)*(sum(v)+1)/2)/(sum(v)*sum(1-v))
stopifnot(abs(roc_auc(h$cases,h$controls)-e)<1e-12,is.na(roc_auc(0,3)))
# Deliberate bin collisions: exact AUC must lie within the ambiguity bound.
h<-roc_hist(roc_bin(c(.100001,.100009)),c(0,1));stopifnot(roc_auc(h$cases,h$controls)==.5)
cat('PASS: count-pooled ROC/AUC matches centralized ranks; ties, patient weights, bin boundaries, empty cases and quantization\n')
