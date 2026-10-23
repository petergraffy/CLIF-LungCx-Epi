source("diagnostics/federated/model_lib.R")
check<-function(x,label){if(!isTRUE(x))stop(label);cat("PASS:",label,"\n")}
for(n in c(10,25,101,500))for(count in 0:n){p<-count/n;b<-fed_interval(p,n);check_n<-sum((0:n)/n>=b[1]-1e-9 & (0:n)/n<=b[2]+1e-9);stopifnot(b[1]<=p+1e-9,b[2]>=p-1e-9,check_n>=10)}
check(TRUE,"fingerprint intervals cover truth and at least ten possible patient counts")
q<-fed_project(c(0,.4,0),c(.4,1,.4));check(abs(sum(q)-1)<1e-8 && all(q>=c(0,.4,0)) && all(q<=c(.4,1,.4)),"interval reconstruction is a feasible probability distribution")
p<-tidyr::expand_grid(hospitalization_id=c("syntheticA","syntheticB"),h=0:71) |>
 dplyr::mutate(support=ifelse(hospitalization_id=="syntheticA","ROOM_AIR","IMV"),hypox="H0",hyper="CU",joint=paste(support,hypox,sep="|"),joint=paste(support,paste(hypox,hyper,sep="_"),sep="|"))
x<-fed_patient_features(p,c("syntheticA","syntheticB"));proto<-x;rownames(proto)<-c("P1","P2")
a<-fed_assign(x,proto,c(72,72));check(identical(a$phenotype,c("P1","P2")),"frozen profile assignment recognizes distinct support trajectories")
a<-fed_assign(x,proto,c(0,72));check(a$phenotype[1]=="I" && a$assignment_status[1]=="insufficient_support","unknown support stays retained under explicit insufficient-information label")
check(inherits(tryCatch(fed_write(data.frame(patient_id="forbidden"),tempfile()),error=identity),"error"),"aggregate export rejects patient identifier columns")
set.seed(7);xx<-cbind(`(Intercept)`=1,z=rnorm(300));yy<-rbinom(300,1,plogis(xx%*%c(-.2,1)));f<-ridge_logistic(xx,yy,paste0("synthetic",1:300))
check(f$converged && f$beta[2]>.5 && all(diag(f$covariance)>0),"local ridge predictor converges with positive clustered variance")
check(auc_rank(c(0,0,1,1),c(.1,.2,.8,.9))==1,"held-out AUC uses correctly ordered predictions")
root<-fed_root()
for(site in c("UCMC","NU","RUSH")){m<-jsonlite::fromJSON(file.path(root,"coordinator",paste0("loo_",site),"frozen_manifest.json"));check(!site%in%m$training_sites && !any(grepl(paste0("/",site,"/"),m$discovery_inputs)),paste("held-out site absent from discovery inputs:",site))}
cat("Federated privacy, assignment, estimation and held-out input checks passed\n")
