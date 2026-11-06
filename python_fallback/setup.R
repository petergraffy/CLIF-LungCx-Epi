# Restore pinned R dependencies without installing, loading or removing Arrow.
source('python_fallback/environment.R')
if(as.character(getRversion())!='4.4.2')stop('Use R 4.4.2; renv does not install R.')
dir.create(clif_bootstrap,recursive=TRUE,showWarnings=FALSE)
dir.create(clif_library,recursive=TRUE,showWarnings=FALSE)
.libPaths(c(clif_library,clif_bootstrap,.libPaths()))
renv_pinned<-function()clif_package_matches('renv','1.1.5')
if(!renv_pinned()){
 for(url in c('https://cloud.r-project.org/src/contrib/renv_1.1.5.tar.gz','https://cloud.r-project.org/src/contrib/Archive/renv/renv_1.1.5.tar.gz')){
  try(utils::install.packages(url,repos=NULL,type='source',lib=clif_bootstrap),silent=TRUE)
  if(renv_pinned())break
 }
}
if(!renv_pinned())stop('Unable to bootstrap renv 1.1.5.')
lock<-renv::lockfile_read('renv.lock');lock$Packages$arrow<-NULL
matching<-names(lock$Packages)[vapply(names(lock$Packages),function(p)clif_package_matches(p,lock$Packages[[p]]$Version),logical(1))]
if(length(matching))renv::hydrate(project=clif_root,library=clif_library,sources=.libPaths(),packages=matching,prompt=FALSE,report=FALSE)
renv::restore(project=clif_root,library=clif_library,lockfile=lock,prompt=FALSE)
.libPaths(c(clif_library,clif_bootstrap,.Library));clif_fallback_check_environment()
cat('Fallback R environment ready. R Arrow was not required. Keep renv.lock unchanged.\n')
