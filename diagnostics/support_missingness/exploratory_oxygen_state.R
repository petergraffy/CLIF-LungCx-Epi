source("code/lib/load.R")
out <- Sys.getenv("CLIF_SUPPORT_AUDIT_OUTPUT", unset = "output/run_UCMC_support_audit_20261022")
a <- readRDS(file.path(init_run(read_run_config()),"local","derivation.rds"))$analyses$h72_primary
e <- readRDS(file.path(out,"oxygen_evidence_local.rds"))
sq <- a$fit$sequence
panel <- tibble::tibble(hospitalization_id=rep(a$fit$ids,ncol(sq)),
 h=rep(seq_len(ncol(sq))-1L,each=nrow(sq)),joint=as.vector(sq)) |>
 dplyr::left_join(e |> dplyr::select(hospitalization_id,h,positive_flow_observed,positive_flow_carried),
  by=c("hospitalization_id","h"))
use <- which(panel$positive_flow_observed | panel$positive_flow_carried)
parts <- strsplit(panel$joint[use],"|",fixed=TRUE)
stopifnot(all(vapply(parts,function(x)x[1]=="UNKNOWN",logical(1))))
panel$joint[use] <- vapply(parts,function(x)paste("OXYGEN_UNSPECIFIED",x[2],sep="|"),character(1))
fit <- sequence_fit(panel,"joint",8L)
fit$distance <- NULL
saveRDS(list(fit=fit,description="Exploratory oxygen evidence state; unchanged primary physiology; not a primary rerun"),
 file.path(out,"exploratory_oxygen_fit_local.rds"))
readr::write_csv(fit$silhouette,file.path(out,"exploratory_oxygen_silhouette.csv"))
map <- tibble::tibble(hospitalization_id=fit$ids,exploratory_cluster=fit$labels)
d <- a$data |> dplyr::left_join(map,by="hospitalization_id")
statehours <- panel |> dplyr::mutate(support=vapply(strsplit(joint,"|",fixed=TRUE),function(x)x[1],character(1))) |>
 dplyr::group_by(hospitalization_id) |> dplyr::summarise(unknown_hours=sum(support=="UNKNOWN"),
 oxygen_unspecified_hours=sum(support=="OXYGEN_UNSPECIFIED"),.groups="drop")
d <- d |> dplyr::left_join(statehours,by="hospitalization_id")
summary <- d |> dplyr::group_by(exploratory_cluster) |> dplyr::summarise(n=dplyr::n(),
 original_unknown_cluster_n=sum(local_cluster==3),unknown_pct=100*mean(unknown_hours/72),
 oxygen_unspecified_pct=100*mean(oxygen_unspecified_hours/72),.groups="drop")
print(summary,width=Inf);print(fit$silhouette)
cat("Agreement with primary:",adjusted_rand(a$fit$labels,fit$labels[match(a$fit$ids,fit$ids)]),"\n")
readr::write_csv(summary,file.path(out,"exploratory_oxygen_clusters.csv"))
readr::write_csv(d |> dplyr::count(local_cluster,exploratory_cluster,name="n"),
 file.path(out,"exploratory_oxygen_cluster_crosswalk.csv"))
stopifnot(length(fit$ids)==nrow(a$data), nrow(panel)==nrow(a$data)*72L,
 sum(statehours$unknown_hours)+sum(statehours$oxygen_unspecified_hours)==sum(a$data$hours_UNKNOWN))
