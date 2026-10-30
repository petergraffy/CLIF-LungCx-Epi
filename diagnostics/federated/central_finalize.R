source("diagnostics/federated/comparator_lib.R")
suppressPackageStartupMessages(library(ggplot2))
root<-fed_root();report<-file.path(root,"reports");for(z in c("tables","figures"))dir.create(file.path(report,z),recursive=TRUE,showWarnings=FALSE)
access<-character()
read_aggregate<-function(path){p<-normalizePath(path,mustWork=TRUE);stopifnot(startsWith(p,paste0(root,"/exchange/"))||startsWith(p,paste0(root,"/coordinator/")));access<<-c(access,p);readr::read_csv(p,show_col_types=FALSE)}
read_manifest<-function(path){p<-normalizePath(path,mustWork=TRUE);stopifnot(startsWith(p,paste0(root,"/exchange/"))||startsWith(p,paste0(root,"/coordinator/")));access<<-c(access,p);jsonlite::fromJSON(p)}
# Copy allowlisted site aggregate results into the exchange before reading.
for(site in c("UCMC","NU","RUSH"))for(family in c("transport","risks","quality")){
 from<-file.path(root,"sites",site,"export",family);to<-file.path(root,"exchange",family,site);dir.create(to,recursive=TRUE,showWarnings=FALSE)
 expected<-switch(family,transport=c("metrics.csv","coverage.csv","manifest.json","extended_metrics.csv","extended_manifest.json"),risks="risks_28day.csv",quality="audit.json")
 stopifnot(setequal(list.files(from),expected));for(name in expected){src<-file.path(from,name);dest<-file.path(to,name);if(file.exists(dest))stopifnot(unname(tools::md5sum(src))==unname(tools::md5sum(dest))) else stopifnot(file.copy(src,dest,overwrite=FALSE))}
}
metrics<-list();coverage<-list();risks<-list();composition<-list();qa<-list();ns<-list();metadata<-list();site_cif<-list()
for(site in c("UCMC","NU","RUSH")){
 metrics[[site]]<-read_aggregate(file.path(root,"exchange","transport",site,"extended_metrics.csv"))
 coverage[[site]]<-read_aggregate(file.path(root,"exchange","transport",site,"coverage.csv"))
 risks[[site]]<-read_aggregate(file.path(root,"exchange","risks",site,"risks_28day.csv"))
 composition[[site]]<-read_aggregate(file.path(root,"exchange","outcomes","full",site,"assignment_counts.csv"))
 metadata[[site]]<-read_aggregate(file.path(root,"exchange","outcomes","full",site,"model_metadata.csv"))
 site_cif[[site]]<-read_aggregate(file.path(root,"exchange","outcomes","full",site,"cumulative_incidence.csv"))
 qa[[site]]<-read_manifest(file.path(root,"exchange","quality",site,"audit.json"));stopifnot(qa[[site]]$status=="passed")
 ns[[site]]<-tibble::tibble(site=site,n=qa[[site]]$candidate_hospitalizations)
}
met<-dplyr::bind_rows(metrics);cov<-dplyr::bind_rows(coverage);risk<-dplyr::bind_rows(risks);comp<-dplyr::bind_rows(composition);n<-dplyr::bind_rows(ns)
pool<-read_aggregate(file.path(root,"coordinator","full","outcomes","pooled_associations.csv"));site_effects<-read_aggregate(file.path(root,"coordinator","full","outcomes","site_associations.csv"))
pro<-read_aggregate(file.path(root,"coordinator","full","prototypes.csv"))
phenotypes<-rownames(readRDS(file.path(root,"coordinator","full","prototypes.rds")))
state_labels<-c(ROOM_AIR="Room-air predominant",LOW_O2="Low-flow oxygen predominant",HFNC="HFNC predominant",CPAP="CPAP predominant",NIPPV="NIPPV predominant",IMV="IMV predominant",OTHER="Other support predominant",UNKNOWN="Unknown support predominant")
dominant<-pro |> dplyr::filter(domain=="support") |> dplyr::group_by(phenotype,state) |> dplyr::summarise(p=mean(probability),.groups="drop") |> dplyr::group_by(phenotype) |> dplyr::slice_max(p,n=1,with_ties=FALSE) |> dplyr::ungroup()
labels<-setNames(paste0(dominant$phenotype,": ",state_labels[dominant$state]),dominant$phenotype)
labels<-c(labels,I="I: Insufficient support information")
fed_write(dominant |> dplyr::mutate(label=labels[phenotype]),file.path(report,"tables","phenotype_labels.csv"))
held_counts<-vapply(c("UCMC","NU","RUSH"),function(site)jsonlite::fromJSON(file.path(root,"coordinator",paste0("loo_",site),"frozen_manifest.json"))$k,numeric(1))
for(pair in list(list(x=met,file="transport_metrics.csv"),list(x=cov,file="transport_coverage.csv"),list(x=risk,file="risks_28day.csv"),list(x=comp,file="phenotype_composition_bands.csv"),list(x=n,file="cohort_counts.csv"),list(x=pool,file="pooled_associations.csv"),list(x=site_effects,file="site_associations.csv"),list(x=dplyr::bind_rows(metadata),file="model_metadata.csv"),list(x=dplyr::bind_rows(site_cif),file="site_cumulative_incidence.csv")))fed_write(pair$x,file.path(report,"tables",pair$file))
auc<-met |> dplyr::filter(metric=="auc",model%in%c("baseline","plus_phenotype","simple_summaries","summaries_plus_phenotype")) |> dplyr::mutate(model=factor(model,levels=c("baseline","plus_phenotype","simple_summaries","summaries_plus_phenotype"),labels=c("Baseline","Baseline + phenotype","Baseline + simple summaries","Baseline + summaries + phenotype")))
p<-ggplot(auc,aes(model,estimate,color=site))+geom_pointrange(aes(ymin=conf_low,ymax=conf_high),position=position_dodge(width=.45))+
 coord_flip()+scale_y_continuous(breaks=seq(.3,1,.1))+labs(title="Prediction transported to a held-out site",subtitle="28-day post-landmark hospital death or hospice | 200 patient bootstrap replicates",x=NULL,y="Held-out AUC with 95% confidence interval",color="Held-out site",caption="Each predictor and its phenotype profiles were learned at the other two sites.\nSimple summaries include support duration, ARF detection and missingness. Comparator added as an exploratory robustness check.")+theme_minimal(base_size=12)+theme(legend.position="bottom",plot.caption=element_text(hjust=0))
ggsave(file.path(report,"figures","held_out_auc.png"),p,width=11,height=5.5,dpi=180)
delta<-met |> dplyr::filter(metric=="auc",model%in%c("plus_minus_baseline","phenotype_beyond_summaries")) |> dplyr::mutate(comparison=ifelse(model=="plus_minus_baseline","Phenotype beyond baseline","Phenotype beyond simple summaries"))
p<-ggplot(delta,aes(site,estimate,color=site))+geom_hline(yintercept=0,linetype=2,color="grey50")+geom_pointrange(aes(ymin=conf_low,ymax=conf_high))+facet_wrap(~comparison)+labs(title="Incremental discrimination depends on the comparator",x="Held-out site",y="Change in AUC (95% patient-bootstrap interval)",caption="Positive values favor adding phenotype. Uncertainty is conditional on frozen profiles and predictors.")+theme_minimal(base_size=12)+theme(legend.position="none",plot.caption=element_text(hjust=0))
ggsave(file.path(report,"figures","incremental_auc.png"),p,width=10,height=5,dpi=180)
support<-pro |> dplyr::filter(domain=="support") |> dplyr::mutate(state=factor(state,levels=SUPPORT_STATES),phenotype=factor(phenotype,levels=phenotypes,labels=labels[phenotypes]))
colors<-c(ROOM_AIR="#B3CDE3",LOW_O2="#377EB8",HFNC="#4DAF4A",CPAP="#984EA3",NIPPV="#6A3D9A",IMV="#E6550D",OTHER="#D8B365",UNKNOWN="#999999")
p<-ggplot(support,aes(bin*6+3,probability,fill=state))+geom_col(width=6,position=position_stack(reverse=TRUE))+facet_wrap(~phenotype,ncol=2)+scale_fill_manual(values=colors,drop=FALSE,labels=c("Room air","Low-flow oxygen","HFNC","CPAP","NIPPV","IMV","Other","Unknown"))+scale_x_continuous(breaks=c(0,24,48,72),limits=c(0,72),expand=c(0,0))+scale_y_continuous(labels=scales::percent,expand=c(0,0))+labs(title=paste(length(phenotypes),"shared respiratory-support / ARF phenotypes"),subtitle=paste0("Outcome-free consensus; held-out-site discovery k: ",paste(held_counts,collapse=" / ")),x="Hours after first ICU arrival",y="Prototype state probability",fill=NULL,caption="Support channel of joint-state prototypes, reconstructed from interval fingerprints.\nSix-hour modal summaries; equal site weighting. These are aggregate templates, not individual trajectories.")+theme_minimal(base_size=12)+theme(legend.position="bottom",panel.spacing.x=grid::unit(1.5,"lines"),plot.caption=element_text(hjust=0))
ggsave(file.path(report,"figures","shared_phenotype_support_profiles.png"),p,width=11,height=7,dpi=180)
phys<-pro |> dplyr::filter(domain%in%c("hypox","hyper")) |> dplyr::mutate(phenotype=factor(phenotype,levels=phenotypes,labels=labels[phenotypes]),domain=factor(domain,levels=c("hypox","hyper"),labels=c("Hypoxemia evidence","Hypercapnia/acidemia evidence")),state=factor(state,levels=c("H0","H1","HU","C0","C1","CU")))
p<-ggplot(phys,aes(bin*6+3,probability,fill=state))+geom_col(width=6,position=position_stack(reverse=TRUE))+facet_grid(phenotype~domain)+scale_fill_manual(values=c(H0="#A6D96A",H1="#D73027",HU="#BBBBBB",C0="#A6CEE3",C1="#762A83",CU="#777777"),labels=c("Hypoxemia negative","Hypoxemia positive","Hypoxemia unknown","Hypercapnia negative","Hypercapnia positive","Hypercapnia unknown"))+scale_x_continuous(breaks=c(0,24,48,72),expand=c(0,0))+scale_y_continuous(labels=scales::percent,expand=c(0,0))+labs(title="ARF evidence within the shared phenotypes",x="Hours after first ICU arrival",y="Prototype state probability",fill=NULL,caption="Original raw-evidence ARF states; unknown does not imply normal physiology.\nInterval reconstruction and six-hour modal summaries smooth the displayed aggregate templates.")+theme_minimal(base_size=11)+theme(legend.position="bottom",panel.spacing.y=grid::unit(.9,"lines"),plot.caption=element_text(hjust=0))
ggsave(file.path(report,"figures","shared_phenotype_arf_profiles.png"),p,width=11,height=10,dpi=170)
forest<-pool |> dplyr::filter(model=="core",endpoint=="death_or_hospice",status=="ok") |> dplyr::mutate(label=labels[phenotype])
p<-ggplot(forest,aes(label,hazard_ratio))+geom_hline(yintercept=1,linetype=2,color="grey50")+geom_pointrange(aes(ymin=hr_low,ymax=hr_high),color="#2166AC")+coord_flip()+scale_y_log10()+labs(title="Pooled association with post-landmark death or hospice",subtitle="Core-adjusted cause-specific hazard ratios through day 28 | reference: P1",x=NULL,y="Hazard ratio (95% modified Hartung–Knapp interval)",caption="Site estimates; REML random effects. All main-cohort patients retained; information-insufficient cases have a separate baseline-hazard stratum.\nCore adjustment omits measured SOFA because of extensive missingness. Associations are conditional and not causal.")+theme_minimal(base_size=11)+theme(plot.caption=element_text(hjust=0))
ggsave(file.path(report,"figures","pooled_composite_associations.png"),p,width=11,height=5,dpi=180)
r<-risk |> dplyr::filter(endpoint=="death_or_hospice",status=="ok") |> dplyr::mutate(phenotype=factor(phenotype,levels=c(phenotypes,"I"),labels=labels[c(phenotypes,"I")]))
p<-ggplot(r,aes(phenotype,risk,color=site))+geom_pointrange(aes(ymin=conf_low,ymax=conf_high),position=position_dodge(width=.4))+coord_flip()+scale_y_continuous(labels=scales::percent,limits=c(0,1))+labs(title="Observed 28-day hospital death or hospice risk",subtitle="Shared phenotype assignments; site-specific risks",x=NULL,y="Risk with 95% patient-clustered interval",color="Site",caption="Small event partitions are withheld. Insufficient information is a data-availability category, not a clinical phenotype.")+theme_minimal(base_size=11)+theme(legend.position="bottom",plot.caption=element_text(hjust=0))
ggsave(file.path(report,"figures","phenotype_absolute_risks.png"),p,width=11,height=6,dpi=180)
# Structural export and train/test separation audit. This is logical process
# isolation on one machine; it does not establish differential privacy.
files<-list.files(file.path(root,"exchange"),recursive=TRUE,full.names=TRUE);stopifnot(all(grepl("\\.(csv|json)$",files)))
for(f in files[grepl("csv$",files)]){z<-read_aggregate(f);stopifnot(!any(names(z)%in%c("patient_id","hospitalization_id","t0","event_time","admission_dttm","discharge_dttm","recorded_dttm","sequence")))}
for(site in c("UCMC","NU","RUSH")){
 version<-paste0("loo_",site);m<-read_manifest(file.path(root,"coordinator",version,"frozen_manifest.json"));p<-read_manifest(file.path(root,"coordinator",version,"outcomes","frozen_predictor.json"));c<-read_manifest(file.path(root,"coordinator",version,"outcomes","frozen_comparator_predictor.json"))
 stopifnot(!site%in%m$training_sites,!site%in%p$training_sites,!site%in%c$training_sites,!any(grepl(paste0("/",site,"/"),c(m$discovery_inputs,p$input_files,c$input_files))))
}
stopifnot(sum(n$n)==3310,all(met$bootstrap_replicates==200))
fed_json(list(status="passed",cohort_n=sum(n$n),patient_identifiers_exported=0,patient_sequences_exported=0,
 all_site_cohort_checks_passed=TRUE,held_out_training_exclusion_passed=TRUE,all_transport_predictions_recomputed_locally=TRUE,
 logical_isolation_only=TRUE,differential_privacy_claimed=FALSE,external_disclosure_review="still required before real between-institution exchange",
 coordinator_report_reads=unique(access)),file.path(report,"federation_audit.json"))
message("Central report tables, six figures and aggregate-only federation audit complete")
