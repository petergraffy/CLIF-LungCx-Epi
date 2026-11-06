# Reconstruct the local feature cache and apply existing prototypes, never learn clusters.
a<-commandArgs(trailingOnly=TRUE);site<-a[1];variant<-a[2]
if(variant=='primary')source('diagnostics/federated/contract.R') else if(variant=='absorbing')source('diagnostics/revision20261022/absorbing_contract.R') else source('diagnostics/revision20261022/variant_contract.R')
source('diagnostics/stage2/rebuild_lib.R')
input<-Sys.getenv(paste0('CLIF_FEDERATED_INPUT_',site));a<-if(variant=='primary')readRDS(file.path(input,'local/waterfall_derivation.rds'))$analyses$h72_waterfall_support_with_raw_arf else revision_input(site)
window<-if(variant=='h24')24 else if(variant=='h48')48 else 72
proto<-readRDS(file.path(fed_root(),'coordinator/full/prototypes.rds'));fm<-jsonlite::fromJSON(file.path(fed_root(),'coordinator/full/frozen_manifest.json'))
x<-fed_patient_features(a$panel,a$data$hospitalization_id);assignment<-fed_assign(x,proto,window-a$data$hours_UNKNOWN)
a$data$local_cluster<-if(variant=='absorbing')sub('^P','D',assignment$phenotype) else assignment$phenotype
fit<-s2_sequence(a$panel,if(variant=='support_only')'support' else 'joint');fit$labels<-a$data$local_cluster[match(fit$ids,a$data$hospitalization_id)];fit$k<-nrow(proto)
home<-file.path(fed_root(),'sites',site,'local');dir.create(home,recursive=TRUE,showWarnings=FALSE)
saveRDS(list(data=a$data,panel=a$panel,fit=fit,input_root=input,consensus_id=fm$consensus_id,grouping='frozen_template',clustering_performed=FALSE),file.path(home,'discovery.rds'))
cat(site,variant,'local data rebuilt and classified against unchanged frozen prototypes; n=',nrow(a$data),'\n')
