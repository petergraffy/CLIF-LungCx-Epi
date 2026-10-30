source("diagnostics/revision20261022/variant_contract.R")
args<-commandArgs(trailingOnly=TRUE);site<-args[1];root<-fed_root();home<-file.path(root,"sites",site)
local<-file.path(home,"local");out<-file.path(home,"export","discovery");dir.create(local,recursive=TRUE,showWarnings=FALSE)
stopifnot(site %in% c("UCMC","NU","RUSH"));cache<-file.path(local,"discovery.rds")
if(file.exists(cache))stop("Fresh rehearsal required; discovery already exists")
input_root<-Sys.getenv(paste0("CLIF_FEDERATED_INPUT_",site),unset=paste0("output/run_",site,"_waterfall_20261022"))
manifest<-jsonlite::fromJSON(file.path(input_root,"local","input_manifest.json"))
a<-revision_input(site)
stopifnot(manifest$site==site,nrow(a$panel)==REV_WINDOW*nrow(a$data),all(a$data$time_from_landmark>0),!any(a$panel$terminal))
d<-a$data;panel<-a$panel
stopifnot(all(d$admission_dttm>=utc("2018-01-01")),all(d$admission_dttm<utc("2025-01-01")),!anyDuplicated(d$hospitalization_id))
wide<-panel |> dplyr::select(hospitalization_id,h,joint) |> tidyr::pivot_wider(names_from=h,values_from=joint,names_prefix="h") |> dplyr::arrange(hospitalization_id)
alphabet<-fed_domains()$joint
sq<-TraMineR::seqdef(as.data.frame(wide[paste0("h",0:(REV_WINDOW-1))]),alphabet=alphabet,right=NA,left=NA,gaps=NA)
sm<-matrix(2,length(alphabet),length(alphabet),dimnames=list(alphabet,alphabet));diag(sm)<-0
message(site," fitting common-cost OM on ",nrow(d)," sequences")
dm<-as.matrix(TraMineR::seqdist(sq,method="OM",indel=1,sm=sm,with.missing=FALSE))
hc<-hclust(as.dist(dm),method="ward.D2");candidates<-2:min(8,nrow(wide)-1)
sil<-tibble::tibble(k=candidates,silhouette=vapply(candidates,function(k)mean(cluster::silhouette(cutree(hc,k),as.dist(dm))[,"sil_width"]),numeric(1)))
k<-sil$k[which.max(sil$silhouette)];labels<-as.integer(cutree(hc,k));ids<-wide$hospitalization_id
fit<-list(ids=ids,labels=labels,distance=dm,k=k,sequence=as.matrix(wide[paste0("h",0:(REV_WINDOW-1))]),field="joint",silhouette=sil)
st<-bootstrap_stability(fit,d$patient_id[match(ids,d$hospitalization_id)],100,20261022L) |>
 dplyr::group_by(local_cluster) |> dplyr::summarise(mean_jaccard=mean(jaccard),.groups="drop")
d$local_cluster<-labels[match(d$hospitalization_id,ids)]
m<-fed_modal(panel) |> dplyr::left_join(d |> dplyr::select(hospitalization_id,local_cluster),by="hospitalization_id")
features<-fed_features();sig<-list();meta<-list()
for(cl in sort(unique(labels))) {
 z<-d[d$local_cluster==cl,];n<-nrow(z);band<-fed_band_n(n);unknown<-mean(z$hours_UNKNOWN/REV_WINDOW);ub<-fed_interval(unknown,n)
 jb<-floor(st$mean_jaccard[st$local_cluster==cl]*10)/10
 meta[[as.character(cl)]]<-tibble::tibble(site=site,local_cluster=cl,n_lower=band[1],n_upper=band[2],
 unknown_lower=ub[1],unknown_upper=ub[2],jaccard_lower=jb,jaccard_upper=pmin(1,jb+.1),small_cluster=n<10)
 counts<-m[m$local_cluster==cl,] |> dplyr::count(domain,bin,state,name="count")
 f<-features |> dplyr::left_join(counts,by=c("domain","bin","state"));f$count[is.na(f$count)]<-0
 intervals<-t(vapply(f$count/n,function(p)fed_interval(p,n),numeric(2)))
 sig[[as.character(cl)]]<-f |> dplyr::select(domain,bin,state,feature,count) |> dplyr::mutate(site=site,local_cluster=cl,denominator=n,lower=intervals[,1],upper=intervals[,2])
}
fed_write(dplyr::bind_rows(sig),file.path(out,"fingerprints.csv"));fed_write(dplyr::bind_rows(meta),file.path(out,"clusters.csv"))
fed_write(sil |> dplyr::mutate(site=site),file.path(out,"silhouette.csv"))
method_files<-c("diagnostics/revision20261022/variant_contract.R","diagnostics/revision20261022/site_discover.R","diagnostics/revision20261022/lib.R")
fed_json(list(schema=FED_SCHEMA,site=site,specification=fed_spec(),spec_hash=fed_hash(fed_spec()),
 method_hash=fed_hash(as.list(tools::md5sum(method_files))),input_processing=manifest$boundary_rule,
 waterfall_source_sha256=jsonlite::fromJSON(file.path(input_root,"local","waterfall_manifest.json"))$waterfall_source_sha256,
 release="interval templates with exact aggregate counts; no patient rows; not differential privacy; local rehearsal"),file.path(out,"manifest.json"))
fit$distance<-NULL
saveRDS(list(data=d,panel=panel,fit=fit,stability=st,input_root=input_root,legacy_adjusted_rand=adjusted_rand(a$fit$labels,labels[match(a$fit$ids,ids)])),cache)
fed_json(list(site=site,n=nrow(d),k=k,all_assigned=TRUE,support_exclusions=0,source=manifest$baseline_root),file.path(local,"discovery_audit.json"))
message(site," local discovery complete: k=",k,"; exported interval fingerprints only")
