source('scripts/deployment_lib.R')
cfg<-deploy_config(Sys.getenv('CLIF_RUN_CONFIG'));root<-deploy_root(cfg)
a<-readRDS(file.path(root,'primary/sites',cfg$site_name,'local/discovery.rds'))
p<-readRDS(file.path(init_run(cfg),'local/prepared.rds'))
# Provisional local-cluster descriptors use the adopted discovery, not raw-baseline labels.
capture_descriptive_tables(p,a$data,cfg,file.path(root,'stage1_descriptors'),a$panel)
capture_site_supplement(cfg,p$root,Sys.getenv('CLIF_RUN_CONFIG'))
