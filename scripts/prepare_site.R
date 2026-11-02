source('code/lib/load.R')
cfg<-read_run_config(Sys.getenv('CLIF_RUN_CONFIG'))
p<-prepare_revision(cfg)
if(!file.exists(file.path(p$root,'local/derivation.rds')))derive_revision(p,cfg)
