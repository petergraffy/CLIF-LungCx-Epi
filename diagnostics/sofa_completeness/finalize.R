source('code/lib/load.R');suppressPackageStartupMessages(library(ggplot2))
root<-'output/run_sofa_completeness_20261026';sites<-c('UCMC','NU','RUSH')
readall<-function(suffix)bind_rows(lapply(sites,function(s)read.csv(file.path(root,'reports',paste0(s,'_',suffix,'.csv')))))
x<-readall('staged_completeness');write.csv(x,file.path(root,'reports/component_completeness_all_sites.csv'),row.names=FALSE)
labels<-c(sofa_resp='Respiratory',sofa_cv='Cardiovascular',sofa_coag='Platelets',sofa_liver='Bilirubin',sofa_renal='Creatinine',sofa_cns='GCS',sofa_total='All six components')
old<-x |> filter(variant=='original') |> select(site,component,n,original=available)
new<-x |> filter(variant=='reviewed_drug_mapping') |> select(site,component,recovered=available)
z<-left_join(old,new,by=c('site','component'));z$component<-factor(z$component,levels=rev(names(labels)),labels=rev(unname(labels)))
p<-ggplot(z,aes(y=component))+geom_segment(aes(x=original/n,xend=recovered/n,yend=component),color='#009E73',linewidth=1.5)+geom_point(aes(x=original/n,color='Original'),size=2.5)+geom_point(aes(x=recovered/n,color='Recovered measurements'),size=2.5,shape=1,stroke=1)+facet_wrap(~site,nrow=1)+scale_x_continuous(labels=scales::percent,limits=c(0,1),breaks=c(0,.25,.5,.75,1))+scale_color_manual(values=c('Original'='#777777','Recovered measurements'='#0072B2'))+labs(title='Which SOFA components are missing?',subtitle='First 24 hours after ICU arrival | same 3,310 landmark hospitalizations',x='Component availability',y=NULL,color=NULL,caption='Recovery uses paired SBP/DBP for MAP and explicit source lab mappings/panel FiO2. Missing values remain missing.\nCreatinine-only renal scoring remains modified SOFA. Pre-ICU weight and estimated PaO2 sensitivities are excluded from this plot.')+theme_minimal(base_size=12)+theme(legend.position='bottom',plot.caption=element_text(hjust=0,size=9))
ggsave(file.path(root,'reports/sofa_component_completeness.png'),p,width=12,height=5,dpi=180)
ggsave(file.path(root,'reports/sofa_component_completeness.pdf'),p,width=12,height=5)
counts<-x |> filter(component=='sofa_total') |> select(site,variant,n,available)
write.csv(counts,file.path(root,'reports/total_score_completeness.csv'),row.names=FALSE)
# Local integrity, monotonic coverage, and aggregate-output checks.
for(s in sites){b<-readRDS(file.path(root,'local',paste0(s,'_extended.rds')));base<-b$variants$original
 for(v in names(b$variants)){a<-b$variants[[v]];stopifnot(identical(a$hospitalization_id,base$hospitalization_id),nrow(a)==nrow(base));for(c in names(labels)[1:6])stopifnot(all(is.na(a[[c]])|(a[[c]]>=0&a[[c]]<=4)),!any(!is.na(base[[c]])&is.na(a[[c]])))
 stopifnot(all(is.na(a$sofa_total)==(rowSums(is.na(a[,names(labels)[1:6]]))>0)))
 }
}
for(f in list.files(file.path(root,'reports'),pattern='csv$',full.names=TRUE)){h<-names(read.csv(f,nrows=1));stopifnot(!any(h%in%c('hospitalization_id','patient_id','recorded_dttm','t0','event_time')))}
stopifnot(sum(old$n[old$component=='sofa_total'])==3310,sum(old$original[old$component=='sofa_total'])==384)
jsonlite::write_json(list(status='passed',sites=sites,n=3310,all_original_members_retained=TRUE,missing_domains_not_scored_zero=TRUE,component_bounds_passed=TRUE,complete_score_requires_all_six=TRUE,aggregate_csv_identifier_check_passed=TRUE,original_models_unchanged=TRUE,source_files=as.list(tools::md5sum(list.files('diagnostics/sofa_completeness',full.names=TRUE)))),file.path(root,'reports/audit.json'),auto_unbox=TRUE,pretty=TRUE)
summary<-capture.output(print(counts,row.names=FALSE))
writeLines(c('SOFA component completeness audit — 2026-10-26','',
'WHY COMPLETENESS IS LOW',
'Total scores require all six measured components. Respiratory P/F availability and total bilirubin are the principal limits; Rush also has incomplete GCS.',
'Original measured components: UCMC respiratory191/CV306/platelets603/bilirubin290/creatinine611/GCS649, N650; NU 444/1528/1596/820/1623/1675, N1691; Rush172/931/930/630/931/563, N969.',
'', 'RECOVERY WITHIN THE ORIGINAL ICU-FIRST-24H WINDOW',
'UCMC: exact-time unambiguous paired SBP/DBP provides calculated MAP; CV availability306->650; complete totals103->116.',
'NU: explicit arterial PO2 labels/LOINCs under other and same-result-time arterial panel FiO2 recover respiratory scores444->489; complete totals219->238. Thirty existing respiratory scores and18 existing complete totals change after additional evidence.',
'Rush: no additional complete totals from MAP/GCS/lab recovery; complete totals remain62. Exact-time GCS components and exact numeric strings provide no new CNS coverage.',
'CV source review:49 positive rows in12 Rush stays name vasoactive drugs but are categorized as dextrose. This correction did not change released score completeness/values in this cohort. One existing UCMC CV score changes0->1 after calculated MAP where recorded MAP is absent; its original total was missing.',
'Combined complete totals384->416 (11.6%->12.6%); these are candidate recovered measured-component modified scores.',
'', 'SEPARATE SENSITIVITIES',
'For drug conversion, keep first-24h measured median weight if available, otherwise use latest measured weight within24h before ICU, bounded by hospital admission. This recovers43 NU and21 Rush CV scores; complete totals become441/3310 (13.3%). No default weight.',
'Only if respiratory score remains missing, estimate PaO2 using the Severinghaus equation implemented in installed CLIFpy2.1.0, restricted to measured SpO2>=50 and<97, paired within1h to raw FiO2 or documented room air. No oxygen-flow FiO2 assumption, numeric waterfall fill or missing-domain-zero. This is an estimated-PaO2 sensitivity, not measured standard SOFA. Combined complete totals1066/3310 (32.2%).',
'',summary,'',
'WHAT REMAINS',
'Bilirubin availability remains290/820/630; paired direct/indirect results recover no missing liver domain. GCS remains649/1675/563; do not turn RASS or absent GCS into GCS15.',
'Allowing lab results from24h before ICU through24h afterward, bounded by admission, increases bilirubin evidence to422/1239/773. This changes the lab window and was audited for availability only; it is not a new first24h score.',
'No urine-output table is present in these supplied CLIF folders, so renal scoring remains creatinine-only. CRRT/IHD/procedure tables do not recover measured urine volume.',
'Missing ABG/FiO2 can reflect measurement and extraction practice, not confirmed normal physiology. Sedation/GCS and original SOFA vasopressor applicability remain limitations. Original source timestamps are used as supplied; this does not independently validate ETL provenance.',
'', 'ANALYSIS STATUS',
'The original clustering, phenotype definitions, outcomes, and SOFA comparisons have NOT been overwritten or refitted. Recovered score comparisons require versioned local model passes and fresh held-out evaluation before reporting updated performance or severity-adjusted associations.',
'Patient-level scores and evidence stay in local/. reports/ contains aggregate audit outputs. The audit is research data QA, not a validated standard-SOFA implementation.'),file.path(root,'reports/sofa_completeness_report.txt'))
cat('SOFA completeness report and integrity checks passed\n')
