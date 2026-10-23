source("diagnostics/federated/comparator_lib.R")
site<-commandArgs(trailingOnly=TRUE)[1];version<-paste0("loo_",site);root<-fed_root();home<-file.path(root,"sites",site)
manifest<-jsonlite::fromJSON(file.path(root,"coordinator",version,"outcomes","frozen_predictor.json"))
cmp<-jsonlite::fromJSON(file.path(root,"coordinator",version,"outcomes","frozen_comparator_predictor.json"))
stopifnot(!site%in%manifest$training_sites,identical(cmp$training_sites,manifest$training_sites),!any(grepl(paste0("/",site,"/"),c(manifest$input_files,cmp$input_files))))
coefficients<-c(manifest$coefficients,cmp$coefficients)
d<-readRDS(file.path(home,"local",version,"outcomes.rds"))$data
phenotypes<-rownames(readRDS(file.path(root,"coordinator",version,"prototypes.rds")))
y<-as.integer(d$time_from_landmark<=28 & d$event_type%in%c("death","hospice"));pred<-list();metrics<-list()
for(kind in names(coefficients)){
 x<-if(kind%in%c("baseline","plus_phenotype"))model_matrix(d,phenotypes,plus=kind=="plus_phenotype") else comparator_matrix(d,phenotypes,plus=kind=="summaries_plus_phenotype")
 b<-unlist(coefficients[[kind]]);stopifnot(setequal(names(b),colnames(x)))
 pred[[kind]]<-plogis(drop(x%*%b[colnames(x)]));metrics[[kind]]<-pred_metrics(y,pred[[kind]])
}
pairs<-list(plus_minus_baseline=c("plus_phenotype","baseline"),simple_minus_baseline=c("simple_summaries","baseline"),phenotype_beyond_summaries=c("summaries_plus_phenotype","simple_summaries"))
set.seed(20261022);patients<-unique(d$patient_id);by<-split(seq_len(nrow(d)),d$patient_id);boot<-list()
for(i in 1:200){ix<-unlist(by[sample(patients,length(patients),replace=TRUE)],use.names=FALSE)
 bm<-lapply(pred,function(p)pred_metrics(y[ix],p[ix]));v<-unlist(lapply(names(bm),function(k)setNames(bm[[k]],paste0(k,"_",names(bm[[k]])))))
 for(pair in names(pairs))for(metric in c("auc","brier"))v[paste(pair,metric,sep="_")]<-bm[[pairs[[pair]][1]]][metric]-bm[[pairs[[pair]][2]]][metric]
 boot[[i]]<-v}
boots<-do.call(rbind,boot);rows<-list()
for(kind in names(metrics))for(metric in names(metrics[[kind]])){
 key<-paste(kind,metric,sep="_");q<-quantile(boots[,key],c(.025,.975),na.rm=TRUE)
 rows[[key]]<-tibble::tibble(site=site,version=version,model=kind,metric=metric,estimate=unname(metrics[[kind]][metric]),conf_low=q[1],conf_high=q[2],bootstrap_replicates=sum(is.finite(boots[,key])))}
for(pair in names(pairs))for(metric in c("auc","brier")){
 key<-paste(pair,metric,sep="_");delta<-metrics[[pairs[[pair]][1]]][metric]-metrics[[pairs[[pair]][2]]][metric];q<-quantile(boots[,key],c(.025,.975),na.rm=TRUE)
 rows[[key]]<-tibble::tibble(site=site,version=version,model=pair,metric=metric,estimate=unname(delta),conf_low=q[1],conf_high=q[2],bootstrap_replicates=sum(is.finite(boots[,key])))}
out<-file.path(home,"export","transport");fed_write(dplyr::bind_rows(rows),file.path(out,"extended_metrics.csv"))
fed_json(list(schema=FED_SCHEMA,site=site,training_sites=manifest$training_sites,held_out_site=site,
 predictor_id=manifest$predictor_id,comparator_predictor_id=cmp$predictor_id,consensus_id=manifest$consensus_id,
 comparator_spec=comparator_spec(),candidate_hospitalizations=nrow(d),all_hospitalizations_evaluated=TRUE,patient_bootstrap_reps=200),file.path(out,"extended_manifest.json"))
saveRDS(list(predictions=pred,y=y,patient_id=d$patient_id,bootstrap=boots),file.path(home,"local",version,"transport_extended_private.rds"))
message(site," extended held-out evaluation complete")
