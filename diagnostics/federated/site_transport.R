source("diagnostics/federated/model_lib.R")
site<-commandArgs(trailingOnly=TRUE)[1];version<-paste0("loo_",site);root<-fed_root();home<-file.path(root,"sites",site)
manifest<-jsonlite::fromJSON(file.path(root,"coordinator",version,"outcomes","frozen_predictor.json"))
stopifnot(!site%in%manifest$training_sites,!any(grepl(paste0("/",site,"/"),manifest$input_files)))
a<-readRDS(file.path(home,"local",version,"outcomes.rds"));d<-a$data
proto<-readRDS(file.path(root,"coordinator",version,"prototypes.rds"));phenotypes<-rownames(proto)
y<-as.integer(d$time_from_landmark<=28 & d$event_type%in%c("death","hospice"));pred<-list();metrics<-list()
for(kind in c("baseline","plus_phenotype")){
 x<-model_matrix(d,phenotypes,plus=kind=="plus_phenotype");b<-unlist(manifest$coefficients[[kind]])
 stopifnot(setequal(names(b),colnames(x)));pred[[kind]]<-plogis(drop(x%*%b[colnames(x)]));metrics[[kind]]<-pred_metrics(y,pred[[kind]])
}
set.seed(20261022);patients<-unique(d$patient_id);by<-split(seq_len(nrow(d)),d$patient_id);boot<-list()
for(i in 1:200){ix<-unlist(by[sample(patients,length(patients),replace=TRUE)],use.names=FALSE)
 b<-pred_metrics(y[ix],pred$baseline[ix]);p<-pred_metrics(y[ix],pred$plus_phenotype[ix]);boot[[i]]<-c(stats::setNames(b,paste0("baseline_",names(b))),stats::setNames(p,paste0("plus_phenotype_",names(p))),delta_auc=p["auc"]-b["auc"],delta_brier=p["brier"]-b["brier"])}
boots<-do.call(rbind,boot);out<-file.path(home,"export","transport");rows<-list()
for(kind in names(metrics))for(metric in names(metrics[[kind]])){
 name<-paste0(kind,"_",metric);q<-quantile(boots[,name],c(.025,.975),na.rm=TRUE)
 rows[[name]]<-tibble::tibble(site=site,version=version,model=kind,metric=metric,estimate=unname(metrics[[kind]][metric]),conf_low=q[1],conf_high=q[2],bootstrap_replicates=sum(is.finite(boots[,name])))}
for(metric in c("auc","brier")){key<-grep(paste0("^delta_",metric),colnames(boots),value=TRUE);delta<-metrics$plus_phenotype[metric]-metrics$baseline[metric];q<-quantile(boots[,key],c(.025,.975),na.rm=TRUE)
 rows[[paste0("delta_",metric)]]<-tibble::tibble(site=site,version=version,model="plus_minus_baseline",metric=metric,estimate=unname(delta),conf_low=q[1],conf_high=q[2],bootstrap_replicates=sum(is.finite(boots[,key])))}
# Post-hoc label comparison is evaluation only, never a training input.
full<-readRDS(file.path(home,"local","full","outcomes.rds"))$data
stopifnot(identical(full$hospitalization_id,d$hospitalization_id))
classified<-d$phenotype!="I";n<-nrow(d);nb<-fed_band_n(n);cb<-fed_band_n(sum(classified))
agreement<-mean(d$phenotype==full$phenotype);ab<-fed_interval(agreement,n)
coverage<-tibble::tibble(site=site,n_lower=nb[1],n_upper=nb[2],classified_lower=cb[1],classified_upper=cb[2],agreement_lower=ab[1],agreement_upper=ab[2])
fed_write(dplyr::bind_rows(rows),file.path(out,"metrics.csv"));fed_write(coverage,file.path(out,"coverage.csv"))
fed_json(list(schema=FED_SCHEMA,site=site,predictor_id=manifest$predictor_id,consensus_id=manifest$consensus_id,
 training_sites=manifest$training_sites,held_out_site=site,patient_bootstrap_reps=200,all_hospitalizations_evaluated=TRUE,
 candidate_hospitalizations=nrow(d),endpoint="28-day in-hospital death/hospice",label_comparison="ordered outcome-free prototype labels; evaluation only"),file.path(out,"manifest.json"))
saveRDS(list(predictions=pred,y=y,patient_id=d$patient_id,bootstrap=boots),file.path(home,"local",version,"transport_private.rds"))
message(site," held-out evaluation complete; no predictions or patient rows exported")
