source('diagnostics/mimic/portable_lib.R')
root<-mimic_root();m<-jsonlite::fromJSON(file.path(root,'coordinator/frozen_predictor.json'));stopifnot(m$spec_hash==fed_hash(portable_spec()),!('MIMIC'%in%m$training_sites),m$coefficient_sha256==digest::digest(file=file.path(root,'coordinator/coefficients.csv'),algo='sha256'),m$prototype_sha256==digest::digest(file=m$prototype_file,algo='sha256'))
b<-readr::read_csv(file.path(root,'coordinator/coefficients.csv'),show_col_types=FALSE);a<-readRDS(file.path(root,'local/derivation.rds'));d<-a$data
rows<-list();pairsout<-list();metadata<-list();private<-list();risks<-list()
for(sub in portable_spec()$subgroups)for(cohort in c('all','measured_sofa'))for(endpoint in c('death_or_hospice','death')){
 z<-if(sub=='coded_lung_cancer')d else d[d$prior_completed_lung_cancer_encounter,]
 if(cohort=='measured_sofa')z<-z[!is.na(z$sofa_total),]
 y<-as.integer(z$time_from_landmark<=28 & if(endpoint=='death')z$event_type=='death' else z$event_type%in%c('death','hospice'))
 nb<-fed_band_n(nrow(z));eb<-fed_band_n(sum(y));key<-paste(sub,cohort,endpoint)
 metadata[[key]]<-tibble::tibble(subgroup=sub,score_cohort=cohort,endpoint=endpoint,n=nrow(z),patients=dplyr::n_distinct(z$patient_id),events_lower=eb[1],events_upper=eb[2],status=if(nrow(z)>=20 && sum(y)>=10 && sum(1-y)>=10)'ok' else 'withheld_small_endpoint_partition')
 if(nrow(z)<20 || sum(y)<10 || sum(1-y)<10)next
 kinds<-if(cohort=='all')c('baseline','phenotype','summaries','summaries_phenotype') else portable_spec()$models
 pred<-list();metrics<-list()
 for(kind in kinds){
  x<-portable_matrix(z,kind);bb<-b[b$cohort==cohort & b$endpoint==endpoint & b$model==kind,];stopifnot(setequal(bb$term,colnames(x)));beta<-setNames(bb$estimate,bb$term)
  pred[[kind]]<-plogis(drop(x%*%beta[colnames(x)]));metrics[[kind]]<-pred_metrics(y,pred[[kind]])
 }
 comparisons<-list(phenotype_beyond_baseline=c('phenotype','baseline'),phenotype_beyond_summaries=c('summaries_phenotype','summaries'))
 if(cohort=='measured_sofa')comparisons$phenotype_beyond_sofa_summaries<-c('sofa_summaries_phenotype','sofa_summaries')
 set.seed(20261023);by<-split(seq_len(nrow(z)),z$patient_id);patients<-names(by);boot<-list()
 for(i in 1:200){ix<-unlist(by[sample(patients,length(patients),replace=TRUE)],use.names=FALSE)
  bm<-lapply(pred,function(p)pred_metrics(y[ix],p[ix]));v<-unlist(lapply(names(bm),function(k)setNames(bm[[k]],paste0(k,'_',names(bm[[k]])))))
  for(pair in names(comparisons))for(metric in c('auc','brier'))v[paste(pair,metric,sep='_')]<-bm[[comparisons[[pair]][1]]][metric]-bm[[comparisons[[pair]][2]]][metric]
  boot[[i]]<-v
 }
 boot<-do.call(rbind,boot)
 for(kind in names(metrics))for(metric in names(metrics[[kind]])){
  name<-paste(kind,metric,sep='_');q<-quantile(boot[,name],c(.025,.975),na.rm=TRUE)
  rows[[paste(key,name)]]<-tibble::tibble(subgroup=sub,score_cohort=cohort,endpoint=endpoint,model=kind,metric=metric,estimate=unname(metrics[[kind]][metric]),conf_low=q[1],conf_high=q[2])
 }
 for(pair in names(comparisons))for(metric in c('auc','brier')){
  name<-paste(pair,metric,sep='_');q<-quantile(boot[,name],c(.025,.975),na.rm=TRUE)
  pairsout[[paste(key,name)]]<-tibble::tibble(subgroup=sub,score_cohort=cohort,endpoint=endpoint,comparison=pair,metric=metric,estimate=unname(metrics[[comparisons[[pair]][1]]][metric]-metrics[[comparisons[[pair]][2]]][metric]),conf_low=q[1],conf_high=q[2])
 }
 private[[key]]<-list(ids=z$hospitalization_id,patient_id=z$patient_id,predictions=pred,y=y,bootstrap=boot)
}
# Separate observed absolute risks, with rare cause partitions masked together.
for(sub in portable_spec()$subgroups){
 z<-if(sub=='coded_lung_cancer')d else d[d$prior_completed_lung_cancer_encounter,]
 for(ph in sort(unique(z$phenotype))){
  zz<-z[z$phenotype==ph,];counts<-c(death=sum(zz$event_type=='death' & zz$time_from_landmark<=28),hospice=sum(zz$event_type=='hospice'&zz$time_from_landmark<=28));safe<-nrow(zz)>=10 && !any(counts>0 & counts<10)
  nb<-fed_band_n(nrow(zz));risks[[paste(sub,ph)]]<-tibble::tibble(subgroup=sub,phenotype=ph,n_lower=nb[1],n_upper=nb[2],death_risk=if(safe)counts['death']/nrow(zz) else NA_real_,hospice_risk=if(safe)counts['hospice']/nrow(zz) else NA_real_,composite_risk=if(sum(counts)>=10 && nrow(zz)-sum(counts)>=10)sum(counts)/nrow(zz) else NA_real_,status=if(safe)'ok' else 'withheld_small_cause_partition')
 }
}
fed_write(dplyr::bind_rows(rows),file.path(root,'aggregate/transport_metrics.csv'));fed_write(dplyr::bind_rows(pairsout),file.path(root,'aggregate/paired_increments.csv'));fed_write(dplyr::bind_rows(metadata),file.path(root,'aggregate/model_denominators.csv'));fed_write(dplyr::bind_rows(risks),file.path(root,'aggregate/absolute_risks.csv'))
saveRDS(private,file.path(root,'local/predictions_private.rds'))
fed_json(list(status='completed',training_sites=m$training_sites,external_site='MIMIC',no_external_refitting=TRUE,cohort_mismatch='POA unavailable; current-coded cohort and prior completed cancer-encounter subgroup',unchanged_phenotypes=TRUE,portable_models='retrained at original three sites without year/POA covariates, not the earlier full-covariate predictor',bootstrap=200,scope='conditional on frozen learner'),file.path(root,'aggregate/evaluation_manifest.json'))
message('MIMIC external prediction complete; no refitting or recalibration')
