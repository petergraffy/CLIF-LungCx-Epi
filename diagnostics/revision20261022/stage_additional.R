source('diagnostics/revision20261022/lib.R');root<-rev_root();groups<-list()
for(s in fed_sites()){
 for(family in c('absorbing','review')){
 from<-if(family=='absorbing')file.path(root,'sensitivities/absorbing_v2/sites',s,'export') else file.path(root,'sites',s,'export')
 to<-file.path(root,'exchange',family,s);dir.create(to,recursive=TRUE,showWarnings=FALSE)
 files<-if(family=='absorbing')c('early_events.csv','absorbing_manifest.json','absorbing_audit.json') else c('severity_source_audit.csv','continuous_denominators.csv','code_source_catalogue.csv','code_source_audit.json')
 for(nm in files){
  target<-file.path(to,nm)
  if(file.exists(target)){stopifnot(unname(tools::md5sum(target))==unname(tools::md5sum(file.path(from,nm))))} else stopifnot(file.copy(file.path(from,nm),target,overwrite=FALSE))
  if(grepl('csv$',nm))groups[[nm]][[s]]<-readr::read_csv(target,show_col_types=FALSE)
 }
 }
}
for(nm in names(groups))fed_write(dplyr::bind_rows(groups[[nm]]),file.path(root,'reports/tables',if(nm=='early_events.csv')'absorbing_early_events.csv' else nm))
