source("code/lib/load.R")
suppressPackageStartupMessages(library(ggplot2))
out<-Sys.getenv("CLIF_WATERFALL_OUTPUT",unset="output/run_UCMC_waterfall_20261022")
r<-readRDS(file.path(out,"local","waterfall_derivation.rds"))
site<-jsonlite::fromJSON(file.path(out,"local","input_manifest.json"))$site
a<-r$analyses$h72_primary
sizes<-a$data |> dplyr::count(local_cluster,name="n")
labels<-stats::setNames(paste0("Cluster ",sizes$local_cluster,"   n = ",sizes$n),sizes$local_cluster)
profile<-a$panel |> dplyr::inner_join(a$data |> dplyr::select(hospitalization_id,local_cluster),by="hospitalization_id") |>
 dplyr::count(local_cluster,h,support,name="n") |>
 tidyr::complete(local_cluster,h,support=SUPPORT_STATES,fill=list(n=0L)) |>
 dplyr::group_by(local_cluster,h) |> dplyr::mutate(proportion=n/sum(n)) |> dplyr::ungroup()
stopifnot(all(abs((profile |> dplyr::group_by(local_cluster,h) |> dplyr::summarise(x=sum(proportion),.groups="drop"))$x-1)<1e-10))
profile$support<-factor(profile$support,levels=SUPPORT_STATES)
colors<-c(ROOM_AIR="#B3CDE3",LOW_O2="#377EB8",HFNC="#4DAF4A",CPAP="#984EA3",NIPPV="#6A3D9A",IMV="#E6550D",OTHER="#D8B365",UNKNOWN="#999999")
p<-ggplot(profile,aes(h+.5,proportion,fill=support))+
 geom_col(width=1,position=position_stack(reverse=TRUE))+
 facet_wrap(~local_cluster,ncol=2,labeller=as_labeller(labels))+
 scale_fill_manual(values=colors,drop=FALSE,labels=c("Room air","Low-flow oxygen","HFNC","CPAP","NIPPV","IMV","Other","Unknown"))+
 scale_x_continuous(breaks=c(0,24,48,72),limits=c(0,72),expand=c(0,0))+
 scale_y_continuous(labels=scales::percent,breaks=c(0,.5,1),expand=c(0,0))+coord_cartesian(ylim=c(0,1))+
 labs(title=paste(site,"joint support and ARF clusters after CLIF waterfall"),
 subtitle=paste0(nrow(a$data)," hospitalizations at the ICU +72-hour landmark | ",a$fit$k," local clusters"),
 x="Hours after first ICU admission",y="Hospitalizations within cluster",fill=NULL,
 caption="Support profiles summarize clustering on joint support and ARF states. Includes inferred and filled states.\nNumeric backfill disabled; processing bounded before each landmark or terminal event. Cluster numbers are specific to this run.")+
 guides(fill=guide_legend(nrow=2,byrow=TRUE))+theme_minimal(base_size=12)+
 theme(legend.position="bottom",panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),
 strip.text=element_text(face="bold"),plot.caption=element_text(hjust=0))
ggsave(file.path(out,"figures","h72_primary_waterfall_support_profiles.png"),p,width=11,height=max(7,2.3*ceiling(a$fit$k/2)+2),dpi=180)
year<-readr::read_csv(file.path(out,"aggregate","unknown_by_year.csv"),show_col_types=FALSE) |>
 dplyr::filter(analysis_id=="h72_primary") |>
 tidyr::pivot_longer(c(original_unknown_hours,waterfall_unknown_hours),names_to="method",values_to="unknown_hours") |>
 dplyr::mutate(unknown_pct=100*unknown_hours/hours,
  method=ifelse(method=="original_unknown_hours","Original six-hour carry","CLIF waterfall"))
p2<-ggplot(year,aes(admission_year,unknown_pct,color=method,group=method))+
 geom_line(linewidth=1)+geom_point(size=2.7)+
 scale_x_continuous(breaks=2018:2024)+scale_y_continuous(limits=c(0,100),labels=function(x)paste0(x,"%"))+
 scale_color_manual(values=c("Original six-hour carry"="#999999","CLIF waterfall"="#2166AC"))+
 labs(title="Waterfall reduces unknown support without changing the cohort",subtitle=paste0(site," 72-hour landmark cohort | ",nrow(a$data)," hospitalizations"),
 x="Hospital admission year",y="Hours with unknown support",color=NULL,
 caption="Waterfall uses documented categories, device/mode inference and forward filling.\nImproved coverage includes imputation and does not establish that the original extract was defective.")+
 theme_minimal(base_size=12)+theme(legend.position="top",panel.grid.minor=element_blank(),plot.caption=element_text(hjust=0))
ggsave(file.path(out,"figures","unknown_support_before_after_waterfall.png"),p2,width=9,height=5.5,dpi=180)
readr::write_csv(profile,file.path(out,"aggregate","h72_primary_support_profile.csv"))
message("Waterfall comparison figures written")
