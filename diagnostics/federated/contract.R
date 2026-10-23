source("code/lib/load.R")
FED_SCHEMA <- "federated-rehearsal-20261022-v1"
fed_spec <- function() list(schema=FED_SCHEMA,window=72L,bin_hours=6L,seed=20261022L,k_max=8L,
 bootstrap_reps=100L,transport_bootstrap=200L,min_cell=10L,minimum_known_support_hours=12L,
 assignment_max_distance=.75,assignment_min_margin=.01,local_om="constant substitution 2; indel 1; Ward.D2",
 exposure="bounded waterfall support with original raw-evidence ARF states",consensus_weight="equal site; size within site",
 fingerprint="interval-valued modal joint/channel distributions and support-change distributions",
 probability_min_width=.05,model_horizon_days=28L,ridge_lambda=1,
 primary_endpoint="28-day in-hospital death or hospice for prediction; post-landmark cause-specific hazards for association")
fed_root <- function() normalizePath(Sys.getenv("CLIF_FEDERATED_OUTPUT",unset="output/run_federated_rehearsal_20261022"),mustWork=FALSE)
fed_hash <- function(x) digest::digest(x,algo="sha256")
fed_write <- function(x,path) {
 forbidden<-intersect(names(x),c("patient_id","hospitalization_id","recorded_dttm","admission_dttm","discharge_dttm","t0","event_time","sequence"))
 if(length(forbidden)||any(vapply(x,inherits,logical(1),"POSIXt")))stop("Forbidden patient-level export")
 dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
 readr::write_csv(x,path,na="")
}
fed_json <- function(x,path) {dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE);jsonlite::write_json(x,path,auto_unbox=TRUE,pretty=TRUE,null="null",digits=16)}
fed_band_n <- function(n) {c(lower=10*floor(n/10),upper=10*floor(n/10)+9)}
fed_interval <- function(p,n) {
 if(n<10)return(c(lower=0,upper=1))
 w<-max(.05,ceiling((10/n)/.05-1e-10)*.05);w<-min(1,w)
 edges<-unique(c(seq(0,1,by=w),1));if(length(edges)>2 && tail(diff(edges),1)<w-1e-9)edges<-edges[-(length(edges)-1)]
 j<-min(findInterval(p,edges,all.inside=TRUE),length(edges)-1)
 c(lower=edges[j],upper=edges[j+1])
}
fed_project <- function(lo,hi) {
 stopifnot(sum(lo)<=1+1e-7,sum(hi)>=1-1e-7)
 mid<-(lo+hi)/2;a<--2;b<-2
 for(i in 1:70){z<-(a+b)/2;if(sum(pmin(hi,pmax(lo,mid-z)))>1)a<-z else b<-z}
 pmin(hi,pmax(lo,mid-(a+b)/2))
}
fed_domains <- function() list(joint=sort(as.vector(outer(SUPPORT_STATES,PHYS_STATES,paste,sep="|"))),
 support=SUPPORT_STATES,hypox=c("H0","H1","HU"),hyper=c("C0","C1","CU"),change=c("down","same","up","unknown"))
fed_mode <- function(x) names(sort(table(x),decreasing=TRUE))[1]
fed_modal <- function(panel) {
 p<-panel |> dplyr::mutate(bin=as.integer(h%/%6)) |>
 dplyr::group_by(hospitalization_id,bin) |> dplyr::summarise(joint=fed_mode(joint),support=fed_mode(support),hypox=fed_mode(hypox),hyper=fed_mode(hyper),.groups="drop")
 ch<-p |> dplyr::arrange(hospitalization_id,bin) |> dplyr::group_by(hospitalization_id) |>
 dplyr::mutate(rank=match(support,SUPPORT_STATES),last=dplyr::lag(rank),change=dplyr::case_when(
 bin==0~NA_character_,support %in% c("UNKNOWN","OTHER")|dplyr::lag(support)%in%c("UNKNOWN","OTHER")~"unknown",
 rank>last~"up",rank<last~"down",TRUE~"same")) |> dplyr::ungroup()
 dplyr::bind_rows(tidyr::pivot_longer(p,c(joint,support,hypox,hyper),names_to="domain",values_to="state"),
 ch |> dplyr::filter(bin>0) |> dplyr::transmute(hospitalization_id,bin,domain="change",state=change))
}
fed_features <- function() {
 dom<-fed_domains();dplyr::bind_rows(lapply(names(dom),function(d)tidyr::expand_grid(domain=d,bin=if(d=="change")1:11 else 0:11,state=dom[[d]]))) |>
 dplyr::mutate(feature=paste(domain,bin,state,sep=":"))
}
fed_weights <- function(features) {
 w<-c(joint=.4,support=.3,hypox=.1,hyper=.1,change=.1)
 sqrt(w[features$domain]/(ifelse(features$domain=="change",11,12)*2))
}
fed_patient_features <- function(panel,ids) {
 f<-fed_features();m<-fed_modal(panel);m$feature<-paste(m$domain,m$bin,m$state,sep=":")
 x<-matrix(0,length(ids),nrow(f),dimnames=list(ids,f$feature))
 x[cbind(match(m$hospitalization_id,ids),match(m$feature,f$feature))]<-1
 x
}
fed_assign <- function(x,proto,known_hours) {
 w<-fed_weights(fed_features());xx<-sweep(x,2,w,"*");pp<-sweep(proto,2,w,"*")
 dist<-outer(rowSums(xx^2),rowSums(pp^2),"+")-2*tcrossprod(xx,pp)
 dist[dist<0]<-0;dist<-sqrt(dist)
 best<-max.col(-dist,ties.method="first");ds<-dist[cbind(seq_len(nrow(x)),best)]
 second<-apply(dist,1,function(z)sort(z,partial=2)[2]);margin<-second-ds
 status<-ifelse(known_hours<fed_spec()$minimum_known_support_hours,"insufficient_support",
 ifelse(ds>fed_spec()$assignment_max_distance,"poor_match",ifelse(margin<fed_spec()$assignment_min_margin,"ambiguous","classified")))
 tibble::tibble(phenotype=ifelse(status=="classified",rownames(proto)[best],"I"),assignment_status=status,distance=ds,margin=margin)
}
