source('code/lib/load.R')
args<-commandArgs(trailingOnly=TRUE)
config_path<-if(length(args))args[1] else 'config/config.json'
cfg<-read_run_config(config_path)
root<-if(length(args)>1)args[2] else supplement_path(cfg)
capture_site_supplement(cfg,root,config_path)
message('Supplement captured in ',file.path(root,'share','supplement'),'; review before returning. Exact diagnostics stay local.')
