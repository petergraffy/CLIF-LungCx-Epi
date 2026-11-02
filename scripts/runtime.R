# Shared site/coordinator environment. Use Rscript --vanilla from the repo root.
clif_root <- normalizePath(getwd(),winslash='/')
if(!file.exists(file.path(clif_root,'renv.lock')))stop('Run from the repository root (renv.lock is missing).')
clif_library <- file.path(clif_root,'renv','library')
clif_bootstrap <- file.path(clif_root,'renv','bootstrap')
.libPaths(c(clif_library,clif_bootstrap,.libPaths()))
clif_check_environment <- function() {
 if(!requireNamespace('jsonlite',quietly=TRUE))stop('Run Rscript --vanilla scripts/setup.R first.')
 lock<-jsonlite::fromJSON('renv.lock',simplifyVector=FALSE)
 if(as.character(getRversion())!=lock$R$Version)stop('Use R ',lock$R$Version,'; current R is ',getRversion(),'. renv restores packages, not R itself.')
 bad<-character()
 for(p in names(lock$Packages)){
  want<-lock$Packages[[p]]$Version
  if(!requireNamespace(p,quietly=TRUE)||utils::packageVersion(p)!=package_version(want))bad<-c(bad,p)
 }
 if(length(bad))stop('R packages missing or different from renv.lock: ',paste(bad,collapse=', '),'. Run scripts/setup.R; do not update the lockfile at a site.')
 invisible(lock)
}
