source('diagnostics/mimic/portable_lib.R')
root<-mimic_root();out<-file.path(root,'coordinator');if(file.exists(file.path(out,'frozen_predictor.json')))stop('Predictor already frozen')
rows<-list();paths<-character()
for(site in portable_spec()$training_sites){
 from<-file.path(root,'training',site,'export');to<-file.path(root,'exchange/training',site);dir.create(to,recursive=TRUE,showWarnings=FALSE)
 for(nm in c('coefficients.csv','metadata.csv','manifest.json')){target<-file.path(to,nm);if(file.exists(target))stopifnot(unname(tools::md5sum(target))==unname(tools::md5sum(file.path(from,nm)))) else stopifnot(file.copy(file.path(from,nm),target,overwrite=FALSE))}
 stopifnot(jsonlite::fromJSON(file.path(to,'manifest.json'))$spec_hash==fed_hash(portable_spec()))
 rows[[site]]<-readr::read_csv(file.path(to,'coefficients.csv'),show_col_types=FALSE);paths<-c(paths,file.path(to,'coefficients.csv'))
}
z<-dplyr::bind_rows(rows) |> dplyr::group_by(cohort,endpoint,model,term) |> dplyr::summarise(k=sum(status=='ok'&is.finite(std_error)&std_error>0),num=sum(ifelse(status=='ok'&is.finite(std_error)&std_error>0,estimate/std_error^2,0),na.rm=TRUE),den=sum(ifelse(status=='ok'&is.finite(std_error)&std_error>0,1/std_error^2,0),na.rm=TRUE),.groups='drop') |> dplyr::mutate(estimate=ifelse(den>0,num/den,0),status=ifelse(k==0,'zero_prior_unestimated',ifelse(k==1,'single_site','multiple_sites'))) |> dplyr::select(-num,-den)
fed_write(z,file.path(out,'coefficients.csv'))
proto<-'output/run_federated_rehearsal_20261022/coordinator/full/prototypes.rds';stopifnot(!any(grepl('/training/MIMIC/',paths,fixed=TRUE)))
fed_json(list(spec=portable_spec(),spec_hash=fed_hash(portable_spec()),training_inputs=paths,training_sites=portable_spec()$training_sites,held_out_site='MIMIC',coefficient_sha256=digest::digest(file=file.path(out,'coefficients.csv'),algo='sha256'),prototype_sha256=digest::digest(file=proto,algo='sha256'),prototype_file=proto,source_md5=as.list(tools::md5sum(c('diagnostics/mimic/portable_lib.R','diagnostics/mimic/site_train_portable.R','diagnostics/mimic/freeze_portable.R')))),file.path(out,'frozen_predictor.json'))
message('Portable predictor frozen using only UCMC/NU/Rush; MIMIC excluded')
