source('diagnostics/apache_adapted/lib.R');suppressPackageStartupMessages(library(ggplot2))
root<-ap_root();report<-file.path(root,'reports');dir.create(file.path(report,'figures'),recursive=TRUE,showWarnings=FALSE);dir.create(file.path(report,'tables'),recursive=TRUE,showWarnings=FALSE)
allowed<-c('component_completeness.csv','component_counts.csv','cohort_counts.csv','provenance_counts.csv','missing_patterns.csv','source_metadata_quality.csv','coefficients.csv','model_metadata.csv','transport_metrics.csv','paired_increments.csv','decision_curves.csv','calibration_bins.csv','evaluation_metadata.csv','component_manifest.json','training_manifest.json','transport_manifest.json','audit.json')
for(site in ap_sites()){
 from<-file.path(root,'sites',site,'export');to<-file.path(root,'exchange',site)
 stopifnot(setequal(list.files(from),allowed));dir.create(to,recursive=TRUE,showWarnings=FALSE)
 for(f in allowed){src<-file.path(from,f);dest<-file.path(to,f);if(file.exists(dest))stopifnot(unname(tools::md5sum(src))==unname(tools::md5sum(dest))) else stopifnot(file.copy(src,dest))}
 audit<-jsonlite::fromJSON(file.path(to,'audit.json'));stopifnot(audit$status=='passed')
}
readall<-function(name)dplyr::bind_rows(lapply(ap_sites(),function(site)readr::read_csv(file.path(root,'exchange',site,name),show_col_types=FALSE)))
for(f in allowed[grepl('[.]csv$',allowed)])ap_export(readall(f),file.path(report,'tables',f))
met<-readall('transport_metrics.csv');inc<-readall('paired_increments.csv');counts<-readall('cohort_counts.csv');components<-readall('component_counts.csv')
labels<-c(complete_hb_proxy='Complete adapted physiology',all_case_missing_indicators='All cases: partial sum + missing indicators',matched_complete_SOFA='Complete for both adapted score and modified SOFA')
models<-c(adapted_only='Adapted score only',baseline_adapted='Baseline + adapted score',adapted_summaries='Baseline + adapted + summaries',adapted_summaries_phenotype='Baseline + adapted + summaries + phenotype',sofa_only='Modified SOFA only',baseline_sofa='Baseline + modified SOFA')
a<-met |> dplyr::filter(endpoint=='death_or_hospice',metric=='auc',model%in%names(models),cohort!='matched_complete_SOFA') |> dplyr::mutate(model=factor(model,levels=names(models),labels=models),cohort=factor(cohort,levels=names(labels),labels=labels))
g<-ggplot(a,aes(model,estimate,color=site))+geom_pointrange(aes(ymin=conf_low,ymax=conf_high),position=position_dodge(width=.55),na.rm=TRUE)+coord_flip()+facet_wrap(~cohort,ncol=1)+scale_y_continuous(limits=c(0,1))+labs(title='Held-out prediction with a CLIF-adapted APACHE-II proxy',subtitle='28-day post-landmark hospital death or hospice',x=NULL,y='AUC (95% patient-bootstrap interval)',color='Held-out site',caption='Hemoglobin-derived hematocrit (3 × Hb); unverified chronic-health points and renal doubling omitted. This is not standard APACHE II or APACHE IV.\nPatients are identical across models within each cohort. Summaries include hour-72 support; bootstrap conditions on frozen other-site learners.')+theme_minimal(base_size=11)+theme(legend.position='bottom',plot.caption=element_text(hjust=0,size=8))
ggsave(file.path(report,'figures/adapted_apache_heldout_auc.png'),g,width=12,height=7,dpi=180)
a<-inc |> dplyr::filter(endpoint=='death_or_hospice',metric=='auc',comparison=='phenotype_beyond_adapted_and_summaries') |> dplyr::mutate(cohort=factor(cohort,levels=names(labels),labels=labels))
g<-ggplot(a,aes(site,estimate,color=cohort))+geom_hline(yintercept=0,linetype=2,color='grey50')+geom_pointrange(aes(ymin=conf_low,ymax=conf_high),position=position_dodge(width=.5),na.rm=TRUE)+coord_flip()+labs(title='Phenotype increment beyond adapted severity and respiratory summaries',x=NULL,y='Change in held-out AUC (95% paired patient-bootstrap interval)',color=NULL,caption='CLIF-adapted APACHE-II APS+age proxy; not standard APACHE II/IV. Zero denotes no increment.\nConditional on frozen original phenotype definitions and other-site coefficient synthesis.')+theme_minimal(base_size=11)+theme(legend.position='bottom',plot.caption=element_text(hjust=0,size=9))+guides(color=guide_legend(nrow=3))
ggsave(file.path(report,'figures/phenotype_increment_beyond_adapted_apache.png'),g,width=12,height=6,dpi=180)
a<-met |> dplyr::filter(endpoint=='death_or_hospice',metric=='auc',cohort=='matched_complete_SOFA',model%in%c('adapted_only','sofa_only','baseline_adapted','baseline_sofa')) |> dplyr::mutate(model=factor(model,levels=names(models),labels=models))
g<-ggplot(a,aes(model,estimate,color=site))+geom_pointrange(aes(ymin=conf_low,ymax=conf_high),position=position_dodge(width=.5),na.rm=TRUE)+coord_flip()+scale_y_continuous(limits=c(0,1))+labs(title='Adapted APACHE proxy and modified SOFA on identical patients',subtitle='Complete for both scores | 28-day post-landmark hospital death or hospice',x=NULL,y='Held-out AUC (95% patient-bootstrap interval)',color='Held-out site',caption='Neither score is a fully observed validated standard score. Matched complete cases are a selected subset.\nPaired AUC/Brier differences, calibration, and decision curves are exported; this plot does not establish universal superiority.')+theme_minimal(base_size=11)+theme(legend.position='bottom',plot.caption=element_text(hjust=0,size=9))
ggsave(file.path(report,'figures/matched_adapted_apache_vs_sofa.png'),g,width=12,height=5,dpi=180)
# Validate original frozen artifact ledger without opening patient data.
ledger<-jsonlite::fromJSON(file.path(ap_primary(),'reports/frozen_artifact_ledger.json'),simplifyVector=FALSE)$ledger
for(x in ledger)stopifnot(unname(tools::md5sum(x$path))==x$md5)
for(site in ap_sites()){
 version<-paste0('loo_',site);m<-jsonlite::fromJSON(file.path(root,'coordinator',version,'frozen_predictor.json'));stopifnot(!site%in%m$training_sites,m$predictor_hash==digest::digest(file=file.path(root,'coordinator',version,'coefficients.csv'),algo='sha256'))
}
files<-list.files(file.path(root,'exchange'),full.names=TRUE,recursive=TRUE)
for(f in files[grepl('csv$',files)]){z<-readr::read_csv(f,show_col_types=FALSE);stopifnot(!any(names(z)%in%c('patient_id','hospitalization_id','t0','event_time','recorded_dttm','sequence')))}
fed_json(list(status='passed',sites=ap_sites(),primary_cohort=sum(counts$n[counts$cohort=='all_case_missing_indicators']),adapted_complete=sum(counts$n[counts$cohort=='complete_hb_proxy']),matched_SOFA_complete=sum(counts$n[counts$cohort=='matched_complete_SOFA']),original_frozen_artifacts_unchanged=TRUE,all_case_sensitivity_retains_every_primary_member=TRUE,held_out_site_exclusion_passed=TRUE,site_predictions_independently_reproduced=TRUE,no_patient_identifiers_in_exchange=TRUE,not_standard_APACHE_II=TRUE,not_APACHE_IV=TRUE),file.path(report,'audit.json'))
fmt<-function(x)if(is.na(x))'unavailable' else sprintf('%.3f',x)
lines<-c('CLIF-adapted APACHE-II comparison sensitivity — October30,2026','',
'Frozen 72h primary phenotypes and original outcomes remain unchanged. New first-ICU-24h adapted APS+age predictions are trained only at the other sites. This is not standard APACHE II or APACHE IV.',
'', 'COHORTS')
for(i in seq_len(nrow(counts)))lines<-c(lines,paste(counts$site[i],counts$cohort[i],counts$n[i],sep=' | '))
lines<-c(lines,'','COMPOSITE HELD-OUT AUC: adapted score only / baseline+adapted / baseline+adapted+summaries / +phenotype')
for(cohort in names(labels))for(site in ap_sites()){
 z<-met[met$cohort==cohort&met$site==site&met$endpoint=='death_or_hospice'&met$metric=='auc',];order<-c('adapted_only','baseline_adapted','adapted_summaries','adapted_summaries_phenotype')
 lines<-c(lines,paste(site,cohort,paste(vapply(z$estimate[match(order,z$model)],fmt,character(1)),collapse=' / '),sep=' | '))
 z<-inc[inc$site==site&inc$cohort==cohort&inc$endpoint=='death_or_hospice'&inc$metric=='auc'&inc$comparison=='phenotype_beyond_adapted_and_summaries',]
 lines<-c(lines,paste0('  phenotype increment: ',fmt(z$estimate),' (',fmt(z$conf_low),', ',fmt(z$conf_high),'); ',z$status))
}
lines<-c(lines,'','MATCHED SAME-PATIENT SCORE COMPARISON: adapted-only AUC / modified-SOFA-only AUC; paired difference')
for(site in ap_sites()){
 z<-met[met$site==site&met$cohort=='matched_complete_SOFA'&met$endpoint=='death_or_hospice'&met$metric=='auc',]
 p<-inc[inc$site==site&inc$cohort=='matched_complete_SOFA'&inc$endpoint=='death_or_hospice'&inc$metric=='auc'&inc$comparison=='adapted_vs_SOFA',]
 lines<-c(lines,paste(site,fmt(z$estimate[z$model=='adapted_only']),fmt(z$estimate[z$model=='sofa_only']),paste0(fmt(p$estimate),' (',fmt(p$conf_low),', ',fmt(p$conf_high),')'),sep=' | '))
}
lines<-c(lines,'','METHODS / LIMITATIONS',
'Independent point-rule implementation corrects units and decimal gaps. All sites use Hb-derived estimated Hct=3*Hb in g/dL for the primary adapted comparator, avoiding site-specific component definitions. Observed-Hct availability scores are descriptive. GCS is valid measured total or exact simultaneous components; sedation status is unverified.',
'Primary adapted totals require all12 observed/substituted APS components and age. Missing domains are unknown. A separate all-case partial-sum plus12 missing-indicator sensitivity is explicitly incomplete, not normal physiology. Chronic-health/surgical context and acute-renal-failure doubling are unverified and excluded; this is APS+age, not full APACHE II.',
'Labs are first24h results after first ICU arrival. Gas collection times are used for specimen identity when valid/in-window, otherwise result time. FiO2 is same-specimen or prior within1h; explicit room air permits0.21. High FiO2 uses sea-level/RQ0.8 A-a approximation. Arterial pH is preferred; absent pH may use serum bicarbonate as a labeled substitution. Invalid units/values are rejected; absent reference units assume canonical CLIF units and are counted.',
'Frozen cached numeric values are retained. NU source moved to CLIF-NU/CLIF_2_1_1: source metadata are attached only on exact hospitalization/time/category/value keys; newer unmatched numeric observations are ignored. Unmatched cached metadata retain explicitly recorded canonical-unit assumptions. Source correspondence counts are exported.',
'Complete and matched-score cohorts are selected subsets. No inference of superiority over standard APACHE, no biological/endotype or causal claims. Existing baseline matrix includes retrospective POA billing features; this is an internal-external development sensitivity, not a validated prospective bedside tool.',
'Coefficientwise inverse-variance ridge synthesis uses a flagged zero prior for unsupported predictive terms. Required score/intercept coefficients must be estimable before evaluation. Performance is released only with>=10 events/non-events; unavailable death estimates remain explicit. AUC/Brier/DCA intervals use200 patient bootstraps conditional on frozen learners; calibration intercept/slope are point estimates. No held-out refitting/recalibration.',
'', 'OUTPUTS','tables/: component/missingness/provenance, same-patient denominators, training coefficients/statuses, held-out AUC/Brier/calibration, paired increments, decision curves, calibration bins.','figures/: held-out adapted score, phenotype increment, same-patient SOFA comparison.','All patient identifiers, component rows, individual predictions, exact denominators and fits stay in site-local local/ directories. Aggregate packages still require site disclosure review.')
writeLines(lines,file.path(report,'adapted_apache_report.txt'))
# Snapshot code and scientific spec for this new sensitivity.
dest<-file.path(root,'provenance');dir.create(dest,recursive=TRUE,showWarnings=FALSE)
sources<-c('code/lib/apache_adapted.R',list.files('diagnostics/apache_adapted',full.names=TRUE),'tests/apache_adapted.R')
for(f in sources){to<-file.path(dest,f);dir.create(dirname(to),recursive=TRUE,showWarnings=FALSE);stopifnot(file.copy(f,to,overwrite=FALSE))}
fed_json(list(specification=ap_spec(),code_md5=as.list(tools::md5sum(sources)),primary_consensus_id=jsonlite::fromJSON(file.path(ap_primary(),'coordinator/full/frozen_manifest.json'))$consensus_id,upstream_helper_not_imported=TRUE),file.path(dest,'manifest.json'))
writeLines(capture.output(sessionInfo()),file.path(dest,'sessionInfo.txt'))
cat('Adapted APACHE report, figures and audits complete\n')
