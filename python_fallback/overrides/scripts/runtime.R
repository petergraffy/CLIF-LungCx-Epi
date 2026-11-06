# Shared site/coordinator environment. Use Rscript --vanilla from the repo root.
clif_root <- normalizePath(getwd(),winslash='/')
if(!file.exists(file.path(clif_root,'renv.lock')))stop('Run from the repository root (renv.lock is missing).')
clif_library <- file.path(Sys.getenv('CLIF_FALLBACK_REPO',unset=clif_root),'renv','library')
clif_bootstrap <- file.path(Sys.getenv('CLIF_FALLBACK_REPO',unset=clif_root),'renv','bootstrap')
.libPaths(c(clif_library,clif_bootstrap,.libPaths()))
# find.package(..., quiet=TRUE) returns character(0) for absent packages.
# Keep DESCRIPTION reads behind an explicit scalar-path check.
clif_package_matches <- function(package, version, finder=find.package) {
 path<-finder(package,quiet=TRUE)
 if(length(path)!=1L||is.na(path)||!nzchar(path))return(FALSE)
 description<-file.path(path,'DESCRIPTION')
 if(!file.exists(description))return(FALSE)
 installed<-tryCatch(read.dcf(description,fields='Version')[1,'Version'],error=function(e)NA_character_)
 if(is.na(installed)||!nzchar(installed))return(FALSE)
 isTRUE(tryCatch(utils::compareVersion(installed,version)==0L,error=function(e)FALSE))
}
clif_check_environment <- function() {
 if(!requireNamespace('jsonlite',quietly=TRUE))stop('Run Rscript --vanilla scripts/setup.R first.')
 lock<-jsonlite::fromJSON('renv.lock',simplifyVector=FALSE)
 if(as.character(getRversion())!=lock$R$Version)stop('Use R ',lock$R$Version,'; current R is ',getRversion(),'. renv restores packages, not R itself.')
 bad<-character()
 for(p in setdiff(names(lock$Packages),"arrow")){
  want<-lock$Packages[[p]]$Version
  if(!requireNamespace(p,quietly=TRUE)||utils::packageVersion(p)!=package_version(want))bad<-c(bad,p)
 }
 if(length(bad))stop('R packages missing or different from renv.lock: ',paste(bad,collapse=', '),'. Run scripts/setup.R; do not update the lockfile at a site.')
 if('arrow'%in%loadedNamespaces())stop('Python fallback must run without loading R Arrow; start a fresh Rscript process.')
 invisible(lock)
}
