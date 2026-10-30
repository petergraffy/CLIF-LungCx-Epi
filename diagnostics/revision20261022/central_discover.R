source("diagnostics/revision20261022/variant_contract.R")
args<-commandArgs(trailingOnly=TRUE);hold<-if(length(args))args[1] else "none"
root<-fed_root();sites<-setdiff(c("UCMC","NU","RUSH"),hold);version<-if(hold=="none")"full" else paste0("loo_",hold)
out<-file.path(root,"coordinator",version);if(file.exists(file.path(out,"frozen_manifest.json")))stop("Consensus already frozen")
features<-fed_features();rows<-list();meta<-list();manifests<-list();read_log<-character()
for(site in sites){dir<-file.path(root,"exchange","discovery",site)
 for(name in c("fingerprints.csv","clusters.csv","manifest.json"))if(!file.exists(file.path(dir,name)))stop("Missing approved aggregate export")
 f<-readr::read_csv(file.path(dir,"fingerprints.csv"),show_col_types=FALSE);m<-readr::read_csv(file.path(dir,"clusters.csv"),show_col_types=FALSE)
 stopifnot(all(c("domain","bin","state","feature","site","local_cluster","lower","upper")%in%names(f)),
 all(f$site==site),all(f$feature %in% features$feature),all(f$lower>=0),all(f$upper<=1),all(f$lower<f$upper))
 manifests[[site]]<-jsonlite::fromJSON(file.path(dir,"manifest.json"));stopifnot(manifests[[site]]$spec_hash==fed_hash(fed_spec()))
 read_log<-c(read_log,file.path(dir,c("fingerprints.csv","clusters.csv","manifest.json")))
 for(cl in sort(unique(f$local_cluster))){z<-f[f$local_cluster==cl,];z<-z[match(features$feature,z$feature),];stopifnot(!anyNA(z$lower))
 p<-numeric(nrow(features))
 for(domain in unique(features$domain))for(bin in unique(features$bin[features$domain==domain])){
 ix<-which(features$domain==domain & features$bin==bin);p[ix]<-fed_project(z$lower[ix],z$upper[ix])}
 key<-paste(site,cl,sep=":");rows[[key]]<-p;meta[[key]]<-m[m$local_cluster==cl,]
 }
}
stopifnot(length(unique(vapply(manifests,function(m)m$method_hash,character(1))))==1)
x<-do.call(rbind,rows);colnames(x)<-features$feature;keys<-dplyr::bind_rows(meta)
keys$key<-names(rows);keys$n_weight<-(keys$n_lower+keys$n_upper)/2
keys$eligible<-!keys$small_cluster & keys$n_lower>=20 & keys$unknown_lower<.5 & keys$jaccard_upper>=.6
train<-which(keys$eligible);if(length(train)<4 || length(unique(keys$site[train]))<2)stop("Insufficient usable cross-site fingerprints")
y<-x[train,,drop=FALSE];ks<-keys[train,]
w<-ks$n_weight/ave(ks$n_weight,ks$site,FUN=sum)/length(unique(ks$site));fw<-fed_weights(features)
distance<-dist(sweep(y,2,fw,"*"));hc<-hclust(distance,method="ward.D2",members=w*length(w))
candidates<-2:min(6,nrow(y)-1);scores<-list();cuts<-list()
for(k in candidates){lab<-as.integer(cutree(hc,k));coverage<-vapply(split(ks$site,lab),function(z)length(unique(z)),integer(1))
 valid<-all(coverage>=2);si<-cluster::silhouette(lab,distance)[,"sil_width"]
 scores[[as.character(k)]]<-tibble::tibble(k=k,weighted_silhouette=sum(w*si),min_contributing_sites=min(coverage),eligible=valid)
 cuts[[as.character(k)]]<-lab}
sc<-dplyr::bind_rows(scores);ok<-which(sc$eligible);if(!length(ok))stop("No candidate phenotype solution represented in >=2 training sites")
best<-ok[which.max(sc$weighted_silhouette[ok])];k<-sc$k[best];lab<-cuts[[as.character(k)]]
proto<-matrix(NA_real_,k,ncol(y))
for(j in seq_len(k)){ix<-which(lab==j);ww<-ks$n_weight[ix]/ave(ks$n_weight[ix],ks$site[ix],FUN=sum)/length(unique(ks$site[ix]));proto[j,]<-colSums(y[ix,,drop=FALSE]*ww)}
colnames(proto)<-features$feature
support_ix<-which(features$domain=="support");rank<-match(features$state[support_ix],SUPPORT_STATES)-1;rank[features$state[support_ix]%in%c("OTHER","UNKNOWN")]<-NA
burden<-apply(proto[,support_ix,drop=FALSE],1,function(z)sum(z*rank,na.rm=TRUE)/sum(z[!is.na(rank)]))
ord<-order(burden);proto<-proto[ord,,drop=FALSE];rownames(proto)<-paste0("P",seq_len(k));lab<-match(lab,ord)
mapdist<-as.matrix(dist(rbind(sweep(x,2,fw,"*"),sweep(proto,2,fw,"*"))))[seq_len(nrow(x)),nrow(x)+seq_len(k),drop=FALSE]
nearest<-max.col(-mapdist,ties.method="first");distmin<-mapdist[cbind(seq_len(nrow(x)),nearest)]
keys$phenotype<-ifelse(keys$small_cluster|keys$unknown_lower>=.5|distmin>.35,"I",rownames(proto)[nearest]);keys$distance<-distmin
keys$used_to_learn_prototypes<-keys$eligible
proto_long<-as.data.frame(t(proto));proto_long$feature<-rownames(proto_long)
proto_long<-tidyr::pivot_longer(proto_long,-feature,names_to="phenotype",values_to="probability") |> dplyr::left_join(features,by="feature")
source_hashes<-lapply(sites,function(s)as.list(tools::md5sum(file.path(root,"exchange","discovery",s,c("fingerprints.csv","clusters.csv","manifest.json")))))
consensus_id<-substr(fed_hash(list(spec=fed_spec(),training_sites=sites,source=source_hashes,prototypes=proto)),1,20)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
fed_write(proto_long,file.path(out,"prototypes.csv"));fed_write(keys,file.path(out,"local_cluster_diagnostic_mapping.csv"));fed_write(sc |> dplyr::mutate(selected=k==.env$k),file.path(out,"candidate_solutions.csv"))
saveRDS(proto,file.path(out,"prototypes.rds"))
fed_json(list(schema=FED_SCHEMA,consensus_id=consensus_id,version=version,training_sites=sites,held_out_site=hold,k=k,
 specification=fed_spec(),spec_hash=fed_hash(fed_spec()),reference="P1",discovery_inputs=read_log,
 source_hashes=source_hashes,assignment="frozen per-hospitalization six-hour modal profile distance; local identifiers stay local",
 coordinator_source=as.list(tools::md5sum("diagnostics/federated/central_discover.R"))),file.path(out,"frozen_manifest.json"))
message(version,": frozen ",k," outcome-free prototypes; training sites ",paste(sites,collapse=","))
print(sc)
