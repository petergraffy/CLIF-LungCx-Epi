source('scripts/runtime.R');source('code/lib/load.R')
base<-tempfile('config_paths_');dir.create(base);base<-normalizePath(base);dir.create(file.path(base,'bin'))
python<-file.path(base,'bin','python')
linked<-suppressWarnings(file.symlink(Sys.which('Rscript'),python))
if(linked){
 cfg<-list(site_name='SYNTHETIC_PATH',repo=base,tables_path=base,python=python,analysis=run_defaults())
 config<-file.path(base,'config.json');jsonlite::write_json(cfg,config,auto_unbox=TRUE,null='null')
 parsed<-read_run_config(config)
 stopifnot(parsed$python==python,parsed$python!=normalizePath(python),file.exists(parsed$python))
 cat('PASS: virtual-environment executable symlinks retain their environment path\n')
}
unlink(base,recursive=TRUE)
