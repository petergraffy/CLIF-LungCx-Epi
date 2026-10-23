source('diagnostics/revision20261022/lib.R')
variant<-Sys.getenv('CLIF_REVISION_VARIANT');stopifnot(nzchar(variant))
REV_WINDOW<-switch(variant,h24=24L,h48=48L,72L)
original_spec<-fed_spec()
fed_spec<-function()modifyList(original_spec,list(schema=paste0('revision-',variant,'-v1'),window=REV_WINDOW,bin_hours=REV_WINDOW/12,minimum_known_support_hours=REV_WINDOW/6,exposure=paste('bounded waterfall support with RAW physiology sensitivity',variant)))
# Twelve bins retain relative order while using 2/4/6-hour resolution.
body(fed_modal)<-parse(text=gsub('h%/%6','h%/%(REV_WINDOW/12)',paste(deparse(body(fed_modal)),collapse='\n'),fixed=TRUE))[[1]]
if(variant=='support_only'){
 fed_domains<-function()list(joint=SUPPORT_STATES,support=SUPPORT_STATES,hypox=c('H0','H1','HU'),hyper=c('C0','C1','CU'),change=c('down','same','up','unknown'))
 fed_weights<-function(features){w<-c(joint=.45,support=.45,hypox=0,hyper=0,change=.1);sqrt(w[features$domain]/(ifelse(features$domain=='change',11,12)*2))}
}
