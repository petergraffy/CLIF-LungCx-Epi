source('diagnostics/revision20261022/score_lib.R')
site<-commandArgs(trailingOnly=TRUE)[1];version<-paste0('loo_',site);root<-rev_root();home<-file.path(root,'sites',site);out<-file.path(home,'export','scores')
b<-readr::read_csv(file.path(root,'coordinator',version,'score_coefficients.csv'),show_col_types=FALSE);m<-jsonlite::fromJSON(file.path(root,'coordinator',version,'score_predictor.json'));stopifnot(!site%in%m$training_sites,m$held_out_site==site,m$spec_hash==fed_hash(score_spec()),identical(digest::digest(file=file.path(root,'coordinator',version,'score_coefficients.csv'),algo='sha256'),m$predictor_hash))
d<-readRDS(file.path(old_root(),'sites',site,'local',version,'outcomes.rds'))$data
rows<-list();decisions<-list();private<-list();paired<-list()
for(cohort in score_spec()$cohorts){
 z<-if(cohort=='measured_sofa')d[!is.na(d$sofa_total),] else d
 for(endpoint in score_spec()$endpoints){
 y<-as.integer(z$time_from_landmark<=28 & if(endpoint=='death')z$event_type=='death' else z$event_type%in%c('death','hospice'));pred<-list();metrics<-list()
 for(kind in score_spec()$models){
 x<-score_matrix(z,kind,cohort);bb<-b[b$cohort==cohort & b$endpoint==endpoint & b$model==kind,];stopifnot(setequal(bb$term,colnames(x)));bb<-setNames(bb$estimate,bb$term)
 pred[[kind]]<-plogis(drop(x%*%bb[colnames(x)]));metrics[[kind]]<-pred_metrics(y,pred[[kind]])
 }
 # Endpoint/event partition must be safe before releasing performance estimates.
 safe<-sum(y)>=10 && sum(1-y)>=10
 pairs<-list(phenotype_beyond_sofa=c('sofa_summaries_phenotype','sofa_summaries'),phenotype_vs_sofa=c('phenotype','sofa_only'),phenotype_beyond_summaries=c('summaries_phenotype','summaries'),phenotype_vs_baseline_sofa=c('phenotype','sofa'))
 set.seed(20261022);by<-split(seq_len(nrow(z)),z$patient_id);patients<-names(by);boot<-list();dcboot<-list();thresholds<-seq(.05,.6,.05)
 for(i in 1:200){ix<-unlist(by[sample(patients,length(patients),replace=TRUE)],use.names=FALSE)
 bm<-lapply(pred,function(p)c(auc=auc_rank(y[ix],p[ix]),brier=mean((y[ix]-p[ix])^2)));v<-unlist(lapply(names(bm),function(k)setNames(bm[[k]],paste0(k,'_',names(bm[[k]])))))
 for(pair in names(pairs))for(metric in c('auc','brier'))v[paste(pair,metric,sep='_')]<-bm[[pairs[[pair]][1]]][metric]-bm[[pairs[[pair]][2]]][metric]
 boot[[i]]<-v
 dcboot[[i]]<-unlist(lapply(pred,function(pr)vapply(thresholds,function(t)mean(y[ix]*(pr[ix]>=t))-mean((1-y[ix])*(pr[ix]>=t))*t/(1-t),numeric(1))))
 }
 boots<-do.call(rbind,boot);db<-do.call(rbind,dcboot)
 for(kind in names(metrics))for(metric in names(metrics[[kind]])){
 key<-paste(kind,metric,sep='_');q<-if(metric%in%c('auc','brier'))quantile(boots[,key],c(.025,.975),na.rm=TRUE) else c(NA,NA)
 rows[[paste(cohort,endpoint,key)]]<-tibble::tibble(site=site,cohort=cohort,endpoint=endpoint,model=kind,metric=metric,estimate=if(safe)unname(metrics[[kind]][metric]) else NA_real_,conf_low=if(safe)q[1] else NA_real_,conf_high=if(safe)q[2] else NA_real_,status=if(safe)'ok' else 'withheld_small_endpoint_partition')
 }
 for(pair in names(pairs))for(metric in c('auc','brier')){
 key<-paste(pair,metric,sep='_');q<-quantile(boots[,key],c(.025,.975),na.rm=TRUE);delta<-metrics[[pairs[[pair]][1]]][metric]-metrics[[pairs[[pair]][2]]][metric]
 paired[[paste(cohort,endpoint,key)]]<-tibble::tibble(site=site,cohort=cohort,endpoint=endpoint,comparison=pair,metric=metric,estimate=if(safe)unname(delta) else NA_real_,conf_low=if(safe)q[1] else NA_real_,conf_high=if(safe)q[2] else NA_real_,status=if(safe)'ok' else 'withheld_small_endpoint_partition')
 }
 for(j in seq_along(pred))for(k in seq_along(thresholds)){
 kind<-names(pred)[j];t<-thresholds[k];v<-pred[[j]]>=t;q<-quantile(db[,(j-1)*length(thresholds)+k],c(.025,.975),na.rm=TRUE)
 decisions[[paste(cohort,endpoint,kind,k)]]<-tibble::tibble(site=site,cohort=cohort,endpoint=endpoint,model=kind,threshold=t,net_benefit=if(safe)mean(y*v)-mean((1-y)*v)*t/(1-t) else NA_real_,conf_low=if(safe)q[1] else NA_real_,conf_high=if(safe)q[2] else NA_real_,treat_all=if(safe)mean(y)-(1-mean(y))*t/(1-t) else NA_real_,treat_none=0)
 }
 private[[paste(cohort,endpoint)]]<-list(ids=z$hospitalization_id,patient_id=z$patient_id,predictions=pred,outcome=y,bootstrap=boots,n=nrow(z),events=sum(y))
 }
}
fed_write(dplyr::bind_rows(rows),file.path(out,'transport_metrics.csv'));fed_write(dplyr::bind_rows(paired),file.path(out,'paired_increments.csv'));fed_write(dplyr::bind_rows(decisions),file.path(out,'decision_curves.csv'));saveRDS(private,file.path(home,'local/score_transport.rds'))
fed_json(list(site=site,held_out_site=site,training_sites=m$training_sites,same_evaluation_patients_across_models=TRUE,held_out_refitting=FALSE,bootstrap=200,bootstrap_scope='conditional on frozen learner; AUC/Brier/DCA intervals; calibration point estimates',score_spec=score_spec()),file.path(out,'transport_manifest.json'))
message(site,' direct modified-SOFA and nested-summary transport complete')
