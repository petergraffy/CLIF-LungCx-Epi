source("code/lib/load.R")
suppressPackageStartupMessages(library(ggplot2))
out <- "output/run_multisite_comparison_20261022"
for(z in c("aggregate","figures","local"))dir.create(file.path(out,z),recursive=TRUE,showWarnings=FALSE)
roots <- c(UCMC="output/run_UCMC_waterfall_20261022",NU="output/run_NU_waterfall_20261022",RUSH="output/run_RUSH_waterfall_20261022")
summaries<-list();years<-list();clusters<-list();profiles<-list();audits<-list();all_analyses<-list();provenance<-list();unknown_sources<-list()
for(site in names(roots)) {
 root<-roots[[site]]; input<-jsonlite::fromJSON(file.path(root,"local","input_manifest.json"))
 b<-readRDS(file.path(input$baseline_root,"local","derivation.rds"))
 p<-readRDS(file.path(input$baseline_root,"local","prepared.rds"))
 w<-readRDS(file.path(root,"local","waterfall_derivation.rds"))
 stopifnot(all(p$cohort$admission_dttm >= as.POSIXct("2018-01-01",tz="UTC")),
 all(p$cohort$admission_dttm < as.POSIXct("2025-01-01",tz="UTC")))
 d<-b$analyses$h72_primary$data
 flow<-landmark_cohort(p$cohort,72) |> dplyr::count(landmark_reason,name="n") |> dplyr::mutate(site=site)
 readr::write_csv(flow,file.path(out,"aggregate",paste0(site,"_landmark_flow.csv")))
 stopifnot(identical(w$baseline_source_id,b$source_id),length(w$analyses)==10)
 wf_manifest<-jsonlite::fromJSON(file.path(root,"local","waterfall_manifest.json"))
 stopifnot(identical(wf_manifest$numeric_bfill,FALSE))
 provenance[[site]]<-list(baseline_root=input$baseline_root,waterfall_root=root,baseline_source_id=b$source_id,waterfall_source_sha256=wf_manifest$waterfall_source_sha256,clifpy_version=wf_manifest$clifpy_version)
 cmp<-w$comparison
 all_analyses[[site]]<-cmp |> dplyr::mutate(site=site)
 raw<-cmp[cmp$analysis_id=="h72_waterfall_support_with_raw_arf",]
 full<-cmp[cmp$analysis_id=="h72_primary",]
 stopifnot(raw$n==nrow(d),full$n==nrow(d),raw$physiology_changed_hours==0)
 summaries[[site]]<-tibble::tibble(site=site,eligible_hospitalizations=nrow(p$cohort),eligible_patients=dplyr::n_distinct(p$cohort$patient_id),
 primary_n=nrow(d),primary_patients=dplyr::n_distinct(d$patient_id),first_year=min(d$admission_year),last_year=max(d$admission_year),
 original_unknown_pct=raw$original_unknown_pct,waterfall_unknown_pct=raw$waterfall_unknown_pct,
 entirely_unknown_before=sum(d$hours_UNKNOWN==72),entirely_unknown_after=sum(w$analyses$h72_primary$data$hours_UNKNOWN==72),
 original_k=raw$original_k,support_waterfall_raw_arf_k=raw$waterfall_k,full_waterfall_k=full$waterfall_k,
 raw_arf_n=sum(d$measured_arf),full_waterfall_arf_n=full$positive_arf_hospitalizations,
 deaths=sum(d$event_type=="death"),death_pct=100*mean(d$event_type=="death"),hospice=sum(d$event_type=="hospice"),hospice_pct=100*mean(d$event_type=="hospice"),
 median_age=median(d$age_years,na.rm=TRUE),
 raw_both_arf_domains_unknown_pct=100*mean(w$analyses$h72_primary$panel$old_physiology=="HU_CU"),
 measured_sofa_n=sum(!is.na(d$sofa_total)),admission_code_unknown_n=sum(d$code_at_admission=="unknown"),
 original_silhouette=raw$original_silhouette,support_waterfall_raw_arf_silhouette=raw$waterfall_silhouette,
 full_waterfall_silhouette=full$waterfall_silhouette,adjusted_rand_support_raw_arf=raw$adjusted_rand)
 u<-w$analyses$h72_primary$data |> dplyr::filter(hours_UNKNOWN==72) |> dplyr::select(hospitalization_id,t0)
 raw_rs<-as.data.frame(arrow::read_parquet(file.path(root,"local","raw_respiratory.parquet")))
 raw_rs$recorded_dttm<-utc(raw_rs$recorded_dttm)
 ur<-raw_rs |> dplyr::inner_join(u,by="hospitalization_id") |> dplyr::filter(recorded_dttm>=t0,recorded_dttm<t0+72*3600)
 unknown_sources[[site]]<-tibble::tibble(site=site,entirely_unknown_n=nrow(u),
 no_respiratory_rows_in_72h_n=nrow(u)-dplyr::n_distinct(ur$hospitalization_id),
 any_known_source_category_n=dplyr::n_distinct(ur$hospitalization_id[support_state(ur$device_category)!="UNKNOWN"]),
 any_positive_source_flow_n=dplyr::n_distinct(ur$hospitalization_id[which(!is.na(ur$lpm_set) & ur$lpm_set>0)]),
 any_source_fio2_n=dplyr::n_distinct(ur$hospitalization_id[which(!is.na(ur$fio2_set))]))
 years[[site]]<-readr::read_csv(file.path(root,"aggregate","unknown_by_year.csv"),show_col_types=FALSE) |> dplyr::filter(analysis_id=="h72_primary") |> dplyr::mutate(site=site)
 audits[[site]]<-readr::read_csv(file.path(root,"aggregate","waterfall_source_conflict_checks.csv"),show_col_types=FALSE) |> dplyr::mutate(site=site)
 a<-w$analyses$h72_waterfall_support_with_raw_arf
 prof<-a$panel |> dplyr::inner_join(a$data |> dplyr::select(hospitalization_id,local_cluster),by="hospitalization_id") |>
 dplyr::count(local_cluster,h,support,name="n") |> tidyr::complete(local_cluster,h,support=SUPPORT_STATES,fill=list(n=0L)) |>
 dplyr::group_by(local_cluster,h) |> dplyr::mutate(proportion=n/sum(n)) |> dplyr::ungroup() |> dplyr::mutate(site=site)
 profiles[[site]]<-prof
 phys<-a$panel |> dplyr::inner_join(a$data |> dplyr::select(hospitalization_id,local_cluster),by="hospitalization_id") |>
 dplyr::group_by(local_cluster) |> dplyr::summarise(hypox_positive_hours_pct=100*mean(hypox=="H1"),hypox_unknown_hours_pct=100*mean(hypox=="HU"),hyper_positive_hours_pct=100*mean(hyper=="C1"),hyper_unknown_hours_pct=100*mean(hyper=="CU"),.groups="drop")
 stability<-readr::read_csv(file.path(root,"aggregate","h72_waterfall_support_with_raw_arf_stability.csv"),show_col_types=FALSE)
 dominant<-prof |> dplyr::group_by(local_cluster,support) |> dplyr::summarise(mean_hours_pct=100*mean(proportion),.groups="drop") |> dplyr::group_by(local_cluster) |> dplyr::slice_max(mean_hours_pct,n=1,with_ties=FALSE) |> dplyr::ungroup() |> dplyr::rename(dominant_support=support)
 clusters[[site]]<-a$data |> dplyr::group_by(local_cluster) |> dplyr::summarise(n=dplyr::n(),death_n=sum(event_type=="death"),death_pct=100*mean(event_type=="death"),hospice_n=sum(event_type=="hospice"),arf_n=sum(measured_arf),unknown_hours_pct=100*mean(hours_UNKNOWN/72),.groups="drop") |>
 dplyr::left_join(phys,by="local_cluster") |> dplyr::left_join(dominant,by="local_cluster") |> dplyr::left_join(stability,by="local_cluster") |> dplyr::mutate(site=site)
 rm(b,p,w);gc()
}
s<-dplyr::bind_rows(summaries); y<-dplyr::bind_rows(years);cl<-dplyr::bind_rows(clusters);pr<-dplyr::bind_rows(profiles)
stopifnot(length(unique(vapply(provenance,function(x)x$waterfall_source_sha256,character(1))))==1)
 jsonlite::write_json(list(sites=provenance,study_start="2018-01-01 00:00:00 UTC",study_end_exclusive="2025-01-01 00:00:00 UTC",every_included_admission_within_window=TRUE,labels="local_unmatched",prognostic_validation="not_repeated_with_fold_contained_waterfall"),file.path(out,"local","comparison_manifest.json"),pretty=TRUE,auto_unbox=TRUE)
 readr::write_csv(dplyr::bind_rows(unknown_sources),file.path(out,"aggregate","entirely_unknown_source_audit.csv"))
 readr::write_csv(dplyr::bind_rows(all_analyses),file.path(out,"aggregate","all_analyses_comparison.csv"))
 readr::write_csv(s,file.path(out,"aggregate","site_comparison.csv"));readr::write_csv(cl,file.path(out,"aggregate","support_waterfall_raw_arf_clusters.csv"));readr::write_csv(y,file.path(out,"aggregate","unknown_by_year.csv"));readr::write_csv(dplyr::bind_rows(audits),file.path(out,"aggregate","source_conflict_checks.csv"))
readr::write_csv(pr,file.path(out,"aggregate","support_waterfall_raw_arf_profiles.csv"))
yl<-y |> tidyr::pivot_longer(c(original_unknown_hours,waterfall_unknown_hours),names_to="method",values_to="unknown") |> dplyr::mutate(pct=100*unknown/hours,method=ifelse(method=="original_unknown_hours","Original six-hour carry","CLIF waterfall"))
p<-ggplot(yl,aes(admission_year,pct,color=method))+geom_line(linewidth=.8)+geom_point(size=2)+facet_wrap(~site,nrow=1)+scale_y_continuous(limits=c(0,100))+scale_x_continuous(breaks=2018:2024)+scale_color_manual(values=c("Original six-hour carry"="#999999","CLIF waterfall"="#2166AC"))+labs(title="Unknown respiratory support across sites",subtitle="Same cohort rules and bounded waterfall processing; 72-hour ICU landmark",x="Admission year",y="Unknown support (% of hours)",color=NULL,caption="Filled or inferred support is distinguished from observed documentation. Site/year denominators differ.")+theme_minimal(base_size=11)+theme(legend.position="top",axis.text.x=element_text(angle=45,hjust=1))
ggsave(file.path(out,"figures","unknown_support_by_site_year.png"),p,width=12,height=4.7,dpi=180)
colors<-c(ROOM_AIR="#B3CDE3",LOW_O2="#377EB8",HFNC="#4DAF4A",CPAP="#984EA3",NIPPV="#6A3D9A",IMV="#E6550D",OTHER="#D8B365",UNKNOWN="#999999")
for(site in names(roots)) {
 pp<-pr |> dplyr::filter(.data$site==.env$site);cc<-cl |> dplyr::filter(.data$site==.env$site)
 labels<-stats::setNames(paste0("Cluster ",cc$local_cluster," | n = ",cc$n),cc$local_cluster)
 pp$support<-factor(pp$support,levels=SUPPORT_STATES)
 fig<-ggplot(pp,aes(h+.5,proportion,fill=support))+geom_col(width=1,position=position_stack(reverse=TRUE))+facet_wrap(~local_cluster,ncol=2,labeller=as_labeller(labels))+scale_fill_manual(values=colors,drop=FALSE,labels=c("Room air","Low-flow oxygen","HFNC","CPAP","NIPPV","IMV","Other","Unknown"))+scale_x_continuous(breaks=c(0,24,48,72),expand=c(0,0))+scale_y_continuous(labels=scales::percent,expand=c(0,0))+labs(title=paste(site,"waterfall support with original ARF states"),subtitle="Local state-sequence clusters; outcomes excluded from clustering",x="Hours after first ICU admission",y="Hospitalizations within cluster",fill=NULL,caption="Cluster numbers have no cross-site correspondence. Includes inferred and filled support; numeric backfill disabled.")+theme_minimal(base_size=11)+theme(legend.position="bottom",panel.grid.minor=element_blank())
 ggsave(file.path(out,"figures",paste0(site,"_support_waterfall_raw_arf_clusters.png")),fig,width=10,height=max(7,ceiling(nrow(cc)/2)*2.2+2),dpi=170)
}
print(s,width=Inf);print(cl,width=Inf)
