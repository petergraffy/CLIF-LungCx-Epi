source('diagnostics/apache_adapted/lib.R');all<-list()
for(site in ap_sites()){
 src<-file.path(ap_root(),'sites',site,'export');dest<-file.path(ap_root(),'exchange',site);dir.create(dest,recursive=TRUE,showWarnings=FALSE)
 for(n in c('coefficients.csv','model_metadata.csv','training_manifest.json','component_manifest.json','component_completeness.csv','provenance_counts.csv','missing_patterns.csv'))stopifnot(file.copy(file.path(src,n),file.path(dest,n),overwrite=FALSE))
 m<-jsonlite::fromJSON(file.path(dest,'training_manifest.json'));stopifnot(m$spec_hash==fed_hash(ap_spec()));all[[site]]<-readr::read_csv(file.path(dest,'coefficients.csv'),show_col_types=FALSE)
}
a<-dplyr::bind_rows(all)
for(hold in ap_sites()){
 version<-paste0('loo_',hold);z<-a[a$site!=hold&a$version==version,]
 b<-z |> dplyr::group_by(cohort,endpoint,model,term) |> dplyr::summarise(contributing_sites=sum(status=='ok'&is.finite(std_error)&std_error>0),num=sum(ifelse(status=='ok'&is.finite(std_error)&std_error>0,estimate/std_error^2,0),na.rm=TRUE),den=sum(ifelse(status=='ok'&is.finite(std_error)&std_error>0,1/std_error^2,0),na.rm=TRUE),.groups='drop') |> dplyr::mutate(estimate=ifelse(den>0,num/den,0),status=ifelse(den>0,'estimated','zero_prior_unestimated')) |> dplyr::select(-num,-den)
 dest<-file.path(ap_root(),'coordinator',version);dir.create(dest,recursive=TRUE,showWarnings=FALSE);ap_export(b,file.path(dest,'coefficients.csv'))
 fed_json(list(held_out_site=hold,training_sites=setdiff(ap_sites(),hold),spec_hash=fed_hash(ap_spec()),predictor_hash=digest::digest(file=file.path(dest,'coefficients.csv'),algo='sha256'),input_files=file.path('exchange',setdiff(ap_sites(),hold),'coefficients.csv')),file.path(dest,'frozen_predictor.json'))
}
cat('Other-site adapted-score predictors frozen\n')
