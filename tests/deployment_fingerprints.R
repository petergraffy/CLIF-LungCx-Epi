# Synthetic aggregate fingerprints; no clinical data and no historical outputs.
source('scripts/deployment_lib.R')
a<-commandArgs(trailingOnly=TRUE);v<-a[1];base<-a[2]
Sys.setenv(CLIF_REVISION_VARIANT=v)
source(if(v=='primary')'diagnostics/federated/contract.R' else if(v=='absorbing')'diagnostics/revision20261022/absorbing_contract.R' else 'diagnostics/revision20261022/variant_contract.R')
f<-fed_features();support<-c('ROOM_AIR','LOW_O2','HFNC','IMV')
for(site in c('ALPHA','BETA_2','GAMMA-3')){
 rows<-list()
 for(cl in 1:4){
  z<-f;states<-c(joint=if(v=='support_only')support[cl] else paste0(support[cl],'|HU_CU'),support=support[cl],hypox='HU',hyper='CU',change='same')
  p<-as.numeric(z$state==states[z$domain]);iv<-t(vapply(p,function(p)fed_interval(p,40),numeric(2)))
  z$count<-40*p;z$site<-site;z$local_cluster<-cl;z$denominator<-40;z$lower<-iv[,1];z$upper<-iv[,2];rows[[cl]]<-z
 }
 out<-file.path(base,site,'discovery',v)
 fed_write(dplyr::bind_rows(rows),file.path(out,'fingerprints.csv'))
 fed_write(tibble::tibble(site=site,local_cluster=1:4,n_lower=40,n_upper=40,unknown_lower=0,unknown_upper=.05,jaccard_lower=.9,jaccard_upper=1,small_cluster=FALSE),file.path(out,'clusters.csv'))
 fed_write(tibble::tibble(k=2:5,silhouette=c(.4,.6,.9,.5),site=site),file.path(out,'silhouette.csv'))
 fed_json(list(site=site,spec_hash=fed_hash(fed_spec()),method_hash=paste0('synthetic-',v)),file.path(out,'manifest.json'))
}
