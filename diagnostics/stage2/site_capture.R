source('scripts/deployment_lib.R');source('diagnostics/apache_adapted/lib.R')
a<-commandArgs(trailingOnly=TRUE);if(length(a)!=4L)stop('Usage: site_capture.R SITE BUNDLE CONFIG WORKSPACE')
site<-a[1];bundle<-a[2];cfg<-deploy_config(a[3]);root<-a[4];return<-file.path(root,'stage2_return')
disc<-readRDS(file.path(root,'primary/sites',site,'local/discovery.rds'));im<-jsonlite::fromJSON(file.path(disc$input_root,'local/input_manifest.json'))
p<-readRDS(file.path(im$baseline_root,'local/prepared.rds'));full<-readRDS(file.path(root,'primary/sites',site,'local/full/outcomes.rds'))$data
cfg$analysis$min_cell<-0L
# Recompute supplements under the common Stage2 code, including corrected APACHE units.
supp<-file.path(root,'supplement_capture');Sys.setenv(CLIF_SUPPLEMENT_OUTPUT=supp)
capture_site_supplement(cfg,im$baseline_root,a[3])
deploy_tree(file.path(supp,'share'),file.path(return,'supplement_package'),c('csv','json','txt'))
ap<-readRDS(file.path(root,'apache/sites',site,'local/components.rds'))
ph<-rownames(readRDS(file.path(root,'primary/coordinator/full/prototypes.rds')))
cm<-deploy_read(file.path(bundle,'bundle_manifest.json'))
desc<-file.path(root,'shared_descriptors');capture_descriptive_tables(p,full,cfg,desc,disc$panel,ph,ap,cm$primary_consensus_id)
deploy_tree(file.path(desc,'share/descriptives'),file.path(return,'descriptives'),c('csv','json'))
for(v in DEPLOY_VARIANTS)deploy_tree(file.path(deploy_analysis_root(root,v),'sites',site,'export'),file.path(return,'analyses',v),c('csv','json'))
deploy_tree(file.path(root,'revision/sites',site,'export'),file.path(return,'review'),c('csv','json'))
deploy_tree(file.path(root,'apache/sites',site,'export'),file.path(return,'apache'),c('csv','json'))
deploy_tree(file.path(root,'four_cox/share'),file.path(return,'four_endpoint_cox'),c('csv','json'))
deploy_tree(file.path(root,'prediction'),file.path(return,'prediction'),c('csv','json'))
deploy_tree(file.path(root,'survival_share'),file.path(return,'survival'),c('csv','json'))
deploy_tree(file.path(root,'input_audit'),file.path(return,'input_audit'),c('csv','json'))
# Exact hourly states, standalone Sankey inputs and full/LOO assignment agreement.
panel<-disc$panel |> dplyr::left_join(full |> dplyr::select(hospitalization_id,phenotype),by='hospitalization_id')
for(domain in c('support','hypox','hyper','joint')){
 z<-panel |> dplyr::count(phenotype,h,state=.data[[domain]],name='n') |> dplyr::group_by(phenotype,h) |> dplyr::mutate(denominator=sum(n),proportion=n/denominator) |> dplyr::ungroup()
 fed_write(z |> dplyr::mutate(site=site),file.path(return,'profiles',paste0('hourly_',domain,'.csv')))
}
initial<-panel[panel$h==0,c('hospitalization_id','support')]
sankey<-full |> dplyr::select(hospitalization_id,phenotype) |> dplyr::left_join(initial,by='hospitalization_id') |> dplyr::count(initial_support=support,phenotype,name='n') |> dplyr::mutate(site=site)
fed_write(sankey,file.path(return,'profiles/initial_support_to_phenotype.csv'))
loo<-readRDS(file.path(root,'primary/sites',site,'local',paste0('loo_',site),'outcomes.rds'))$data
# Match the LOO labels to full labels using prototype distance (small k exhaustive permutation).
afull<-readRDS(file.path(root,'primary/coordinator/full/prototypes.rds'));alo<-readRDS(file.path(root,'primary/coordinator',paste0('loo_',site),'prototypes.rds'))
source('diagnostics/federated/contract.R')
if(nrow(afull)==nrow(alo)){
 permutations<-function(v)if(length(v)==1L)matrix(v,1) else do.call(rbind,lapply(seq_along(v),function(i)cbind(v[i],permutations(v[-i]))))
 perms<-permutations(seq_len(nrow(afull)));fw<-fed_weights(fed_features());alo<-alo[,colnames(afull),drop=FALSE]
 cost<-apply(perms,1,function(ix)sum((sweep(alo,2,fw,'*')-sweep(afull[ix,,drop=FALSE],2,fw,'*'))^2))
 ix<-perms[which.min(cost),];mapping<-setNames(rownames(afull)[ix],rownames(alo));mapping<-c(mapping,I='I')
 li<-match(full$hospitalization_id,loo$hospitalization_id)
 agreement<-tibble::tibble(full=full$phenotype,loo_aligned=unname(mapping[loo$phenotype[li]])) |> dplyr::count(full,loo_aligned,name='n') |> dplyr::mutate(site=site)
 fed_write(agreement,file.path(return,'profiles/phenotype_transport_agreement.csv'))
}else fed_write(tibble::tibble(site=site,status='different_k_compare_prototype_distances',full_k=nrow(afull),loo_k=nrow(alo)),file.path(return,'profiles/phenotype_transport_agreement.csv'))
deploy_copy('renv.lock',file.path(return,'environment/renv.lock.json'));writeLines(capture.output(sessionInfo()),file.path(return,'environment/R_sessionInfo.txt'))
