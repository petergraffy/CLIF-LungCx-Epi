source("diagnostics/federated/contract.R")
model_spec<-function()list(horizon_days=28L,core=c("age10","male","black","other_or_unknown_race","year","charlson5","advanced","code_nonfull","code_unknown"),
 cox="cause-specific; censor at 28 days; separate baseline stratum for insufficient-information group; patient-clustered SE",
 severity_sensitivity="add measured-component SOFA; complete cases explicitly reported",prediction="ridge logistic; lambda=1 slopes, unpenalized intercept; 28-day hospital death/hospice",prediction_pool="coefficientwise inverse-variance synthesis; unsupported coefficients retain zero prior and are flagged")
base_matrix<-function(d){race<-tolower(d$race_category);code<-d$code_at_admission
 x<-cbind(`(Intercept)`=1,age10=(d$age_years-70)/10,male=as.integer(tolower(d$sex_category)=="male"),
 black=as.integer(grepl("black|african",race)),other_or_unknown_race=as.integer(!grepl("black|african|white",race)),
 year=d$admission_year-2021,charlson5=d$charlson_noncancer/5,advanced=as.numeric(d$advanced_cancer),
 code_nonfull=as.integer(!code%in%c("full","unknown")),code_unknown=as.integer(code=="unknown"))
 stopifnot(!anyNA(x));x}
model_matrix<-function(d,phenotypes,plus=FALSE){x<-base_matrix(d)
 if(plus){p<-c(phenotypes[-1],"I");z<-vapply(p,function(a)as.numeric(d$phenotype==a),numeric(nrow(d)));colnames(z)<-paste0("phenotype_",p);x<-cbind(x,z)}
 x}
ridge_logistic<-function(x,y,patient_id,lambda=1){
 penalty<-diag(c(0,rep(lambda,ncol(x)-1)));b<-rep(0,ncol(x));b[1]<-qlogis(mean(y))
 objective<-function(beta){eta<-drop(x%*%beta);sum(y*eta-pmax(eta,0)-log1p(exp(-abs(eta))))-.5*sum(beta*(penalty%*%beta))}
 converged<-FALSE
 for(i in 1:100){eta<-drop(x%*%b);p<-plogis(eta);w<-pmax(p*(1-p),1e-8);h<-crossprod(x,x*w)+penalty
 step<-solve(h,drop(crossprod(x,y-p))-penalty%*%b);scale<-1
 while(objective(b+scale*step)<objective(b)-1e-9 && scale>1e-7)scale<-scale/2
 b<-drop(b+scale*step);if(max(abs(scale*step))<1e-7){converged<-TRUE;break}}
 p<-plogis(drop(x%*%b));h<-crossprod(x,x*pmax(p*(1-p),1e-8))+penalty;bread<-solve(h)
 scores<-rowsum(x*(y-p),patient_id,reorder=FALSE);v<-bread%*%crossprod(scores)%*%bread
 names(b)<-colnames(x);dimnames(v)<-list(colnames(x),colnames(x))
 list(beta=b,covariance=v,converged=converged,iterations=i)
}
auc_rank<-function(y,p){n1<-sum(y);n0<-length(y)-n1;if(!n1||!n0)return(NA_real_);(sum(rank(p)[y==1])-n1*(n1+1)/2)/(n1*n0)}
pred_metrics<-function(y,p){p<-pmin(1-1e-8,pmax(1e-8,p));lp<-qlogis(p)
 ci<-tryCatch(coef(suppressWarnings(glm(y~1+offset(lp),family=binomial())))[1],error=function(e)NA_real_)
 cs<-tryCatch(coef(suppressWarnings(glm(y~lp,family=binomial())))[2],error=function(e)NA_real_)
 c(auc=auc_rank(y,p),brier=mean((y-p)^2),calibration_intercept=unname(ci),calibration_slope=unname(cs))}
