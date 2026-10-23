source('diagnostics/revision20261022/score_lib.R')
root<-rev_root();all<-list()
for(site in c('UCMC','NU','RUSH')){
 from<-file.path(root,'sites',site,'export','scores');to<-file.path(root,'exchange','scores',site);dir.create(to,recursive=TRUE,showWarnings=FALSE)
 for(nm in c('coefficients.csv','model_metadata.csv','manifest.json'))stopifnot(file.copy(file.path(from,nm),file.path(to,nm),overwrite=FALSE))
 m<-jsonlite::fromJSON(file.path(to,'manifest.json'));stopifnot(m$spec_hash==fed_hash(score_spec()));all[[site]]<-readr::read_csv(file.path(to,'coefficients.csv'),show_col_types=FALSE)
}
a<-dplyr::bind_rows(all)
for(hold in c('UCMC','NU','RUSH')){
 version<-paste0('loo_',hold);z<-a[a$site!=hold & a$version==version,]
 b<-z |> dplyr::group_by(cohort,endpoint,model,term) |> dplyr::summarise(contributing_sites=sum(status=='ok' & is.finite(std_error)&std_error>0),num=sum(ifelse(status=='ok' & is.finite(std_error)&std_error>0,estimate/std_error^2,0),na.rm=TRUE),den=sum(ifelse(status=='ok'&is.finite(std_error)&std_error>0,1/std_error^2,0),na.rm=TRUE),.groups='drop') |> dplyr::mutate(estimate=ifelse(den>0,num/den,0),status=ifelse(den>0,'estimated','zero_prior_unestimated')) |> dplyr::select(-num,-den)
 out<-file.path(root,'coordinator',version);fed_write(b,file.path(out,'score_coefficients.csv'));fed_json(list(training_sites=setdiff(c('UCMC','NU','RUSH'),hold),held_out_site=hold,predictor_hash=digest::digest(file=file.path(out,'score_coefficients.csv'),algo='sha256'),spec_hash=fed_hash(score_spec()),input_files=file.path(root,'exchange','scores',setdiff(c('UCMC','NU','RUSH'),hold),'coefficients.csv')),file.path(out,'score_predictor.json'))
}
