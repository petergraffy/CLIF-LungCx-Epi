# Restore the pinned project library without changing the global R library.
source('scripts/runtime.R')
if(as.character(getRversion())!='4.4.2')stop('Install/use R 4.4.2 before setup; renv does not install R.')
dir.create(clif_bootstrap,recursive=TRUE,showWarnings=FALSE)
renv_pinned<-function(){p<-find.package('renv',quiet=TRUE);nzchar(p)&&read.dcf(file.path(p,'DESCRIPTION'))[1,'Version']=='1.1.5'}
if(!renv_pinned()){
 urls<-c('https://cloud.r-project.org/src/contrib/renv_1.1.5.tar.gz','https://cloud.r-project.org/src/contrib/Archive/renv/renv_1.1.5.tar.gz')
 for(url in urls){
  try(utils::install.packages(url,repos=NULL,type='source',lib=clif_bootstrap),silent=TRUE)
  if(renv_pinned())break
 }
}
if(!requireNamespace('renv',quietly=TRUE)||as.character(utils::packageVersion('renv'))!='1.1.5')stop('Unable to bootstrap pinned renv 1.1.5; inspect installation/network errors.')
dir.create(clif_library,recursive=TRUE,showWarnings=FALSE)
lock<-renv::lockfile_read(file.path(clif_root,'renv.lock'))
matching<-names(lock$Packages)[vapply(names(lock$Packages),function(p){
 installed<-find.package(p,quiet=TRUE)
 nzchar(installed)&&package_version(read.dcf(file.path(installed,'DESCRIPTION'))[1,'Version'])==package_version(lock$Packages[[p]]$Version)
},logical(1))]
if(length(matching))renv::hydrate(project=clif_root,library=clif_library,sources=.libPaths(),packages=matching,prompt=FALSE,report=FALSE)
renv::restore(project=clif_root,library=clif_library,lockfile=file.path(clif_root,'renv.lock'),prompt=FALSE)
.libPaths(c(clif_library,clif_bootstrap,.Library));clif_check_environment()
cat('Pinned R environment ready. Next: set up Python and run the site preflight.\n')
