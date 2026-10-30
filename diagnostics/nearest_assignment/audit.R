source('diagnostics/federated/contract.R')
root<-'output/run_nearest_assignment_20261030';summaries<-list();transitions<-list();template_checks<-list()
jobs<-c(primary='primary',stats::setNames(paste0('revision/sensitivities/',c('h24','h48','onset','support_only','low_flow','carry0','carry12','raw_support')),c('h24','h48','onset','support_only','low_flow','carry0','carry12','raw_support')))
for(job in names(jobs)){
 source('diagnostics/federated/contract.R')
 if(job!='primary'){Sys.setenv(CLIF_REVISION_VARIANT=job);source('diagnostics/revision20261022/variant_contract.R')}
 current<-file.path(root,jobs[[job]]);parent<-jsonlite::fromJSON(file.path(current,'nearest_assignment_parent.json'))$parent_root
 features<-fed_features();w<-fed_weights(features)
 for(v in c('full','loo_UCMC','loo_NU','loo_RUSH')){
  proto<-readRDS(file.path(current,'coordinator',v,'prototypes.rds'))
  stopifnot(identical(proto,readRDS(file.path(parent,'coordinator',v,'prototypes.rds'))))
  template_checks[[paste(job,v)]]<-TRUE
  pp<-sweep(proto,2,w,'*')
  for(site in c('UCMC','NU','RUSH')){
   d<-readRDS(file.path(current,'sites',site,'local',v,'outcomes.rds'))$data
   old<-readRDS(file.path(parent,'sites',site,'local',v,'outcomes.rds'))$data
   x<-readRDS(file.path(current,'sites',site,'local/patient_features.rds'));xx<-sweep(x,2,w,'*')
   sq<-outer(rowSums(xx^2),rowSums(pp^2),'+')-2*tcrossprod(xx,pp);sq[sq<0]<-0;dm<-sqrt(sq)
   best<-rownames(proto)[max.col(-dm,ties.method='first')]
   sufficient<-fed_spec()$window-d$hours_UNKNOWN>=fed_spec()$minimum_known_support_hours
   stopifnot(identical(old$hospitalization_id,d$hospitalization_id),identical(rownames(x),d$hospitalization_id),identical(old$event_type,d$event_type),identical(old$event_time,d$event_time),
    all(d$phenotype[sufficient]==best[sufficient]),all(d$phenotype[!sufficient]=='I'),all(d$assignment_status[sufficient]=='classified'),all(d$assignment_status[!sufficient]=='insufficient_support'),
    !any(d$assignment_status%in%c('ambiguous','poor_match')),all(d$phenotype[old$phenotype!='I']==old$phenotype[old$phenotype!='I']))
   summaries[[paste(job,v,site)]]<-tibble::tibble(analysis=job,version=v,site=site,n=nrow(d),sufficient_n=sum(sufficient),classified_n=sum(d$phenotype!='I'),insufficient_support_n=sum(!sufficient),ambiguous_n=0L,poor_match_n=0L,newly_classified_n=sum(old$phenotype=='I'&d$phenotype!='I'))
   if(job=='primary'&&v=='full')transitions[[site]]<-tibble::tibble(site=site,previous=old$phenotype,previous_status=old$assignment_status,current=d$phenotype) |> dplyr::count(site,previous,previous_status,current,name='n')
  }
 }
}
for(site in c('UCMC','NU','RUSH')){
 x<-readRDS(file.path(root,'revision/sensitivities/absorbing_v2/sites',site,'local/absorbing_summary.rds'))
 stopifnot(!any(x$assignment$assignment_status%in%c('ambiguous','poor_match')),all(x$assignment$phenotype[x$assignment$assignment_status=='classified']==x$assignment$closest_phenotype[x$assignment$assignment_status=='classified']))
}
dir.create(file.path(root,'reports'),recursive=TRUE,showWarnings=FALSE)
readr::write_csv(dplyr::bind_rows(summaries),file.path(root,'reports/assignment_audit.csv'))
readr::write_csv(dplyr::bind_rows(transitions),file.path(root,'reports/assignment_transitions.csv'))
jsonlite::write_json(list(status='passed',site_version_passes=length(summaries),unchanged_template_versions=length(template_checks),no_ambiguous_or_poor_match_assignments=TRUE,nearest_distances_independently_recomputed=TRUE,cohort_and_raw_outcomes_unchanged=TRUE,absorbing_assignment_passed=TRUE),file.path(root,'reports/assignment_audit.json'),auto_unbox=TRUE,pretty=TRUE)
cat('All108 primary/sensitivity site-version assignments independently match nearest templates\n')
