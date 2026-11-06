# The fallback uses the same lockfile except that R Arrow is not a dependency.
source('scripts/runtime.R')
clif_fallback_check_environment<-function(){
 if(!requireNamespace('jsonlite',quietly=TRUE))stop('Run Rscript --vanilla python_fallback/setup.R first.')
 lock<-jsonlite::fromJSON('renv.lock',simplifyVector=FALSE)
 if(as.character(getRversion())!=lock$R$Version)stop('Use R ',lock$R$Version,'.')
 packages<-setdiff(names(lock$Packages),'arrow')
 bad<-packages[!vapply(packages,function(p)clif_package_matches(p,lock$Packages[[p]]$Version),logical(1))]
 if(length(bad))stop('Restore the fallback R environment first; packages: ',paste(bad,collapse=', '))
 invisible(lock)
}
