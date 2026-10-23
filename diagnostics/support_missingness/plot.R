suppressPackageStartupMessages(library(ggplot2))
out <- Sys.getenv("CLIF_SUPPORT_AUDIT_OUTPUT", unset = "output/run_UCMC_support_audit_20261022")
if(!file.exists(file.path(out,"unknown_by_year_location.csv"))) {
 h <- readRDS(file.path(out,"hour_audit_local.rds"))
 d <- readRDS(file.path(init_run(read_run_config()),"local","derivation.rds"))$analyses$h72_primary$data
 annual <- h |> dplyr::left_join(d |> dplyr::select(hospitalization_id,admission_year),by="hospitalization_id") |>
  dplyr::group_by(admission_year,location) |> dplyr::summarise(hours=dplyr::n(),
   unknown_hours=sum(support=="UNKNOWN"),unknown_pct=100*mean(support=="UNKNOWN"),.groups="drop")
 readr::write_csv(annual,file.path(out,"unknown_by_year_location.csv"))
}
z <- readr::read_csv(file.path(out,"unknown_by_year_location.csv"),show_col_types=FALSE)
z <- z[z$location %in% c("icu","ward"),]
z$location <- factor(z$location,levels=c("icu","ward"),labels=c("ICU","Ward"))
p <- ggplot(z,aes(admission_year,unknown_pct,color=location,group=location))+
 geom_line(linewidth=1)+geom_point(size=2.8)+
 geom_text(data=z[z$admission_year %in% c(2018,2020,2024),],
  aes(label=paste0(round(unknown_pct,1),"%")),vjust=-.9,size=3.5,show.legend=FALSE)+
 scale_color_manual(values=c(ICU="#2166AC",Ward="#B35806"))+
 scale_x_continuous(breaks=2018:2024)+scale_y_continuous(limits=c(0,100),breaks=seq(0,100,20),labels=function(x)paste0(x,"%"))+
 labs(title="Unknown respiratory support rises sharply after 2019",
 subtitle="UCMC 72-hour landmark cohort: 650 hospitalizations; 46,800 hourly states",
 x="Hospital admission year",y="Hours with unknown support",color=NULL,
 caption="Unknown support is retained, not labeled room air. ICU location takes precedence in hours spanning a transfer.\nSupport uses recorded device categories with at most six hours of carry-forward.")+
 theme_minimal(base_size=12)+theme(legend.position="top",panel.grid.minor=element_blank(),
 plot.title=element_text(face="bold"),plot.caption=element_text(hjust=0))
ggsave(file.path(out,"unknown_support_by_year_location.png"),p,width=9,height=5.4,dpi=180)
