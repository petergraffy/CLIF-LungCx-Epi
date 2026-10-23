source("diagnostics/federated/comparator_lib.R")
root<-fed_root()
for(version in c("full","loo_UCMC","loo_NU","loo_RUSH")){
 if(file.exists(file.path(root,"coordinator",version,"outcomes","frozen_comparator_predictor.json")))stop("Comparator predictor already frozen")
 m<-jsonlite::fromJSON(file.path(root,"coordinator",version,"frozen_manifest.json"));rows<-list();paths<-character()
 for(site in m$training_sites){from<-file.path(root,"sites",site,"export","comparators",version);to<-file.path(root,"exchange","comparators",version,site)
 dir.create(to,recursive=TRUE,showWarnings=FALSE)
 for(name in c("coefficients.csv","manifest.json"))stopifnot(file.copy(file.path(from,name),file.path(to,name),overwrite=FALSE))
 # Read only the staged export after copy; no private site files are accessed.
 mm<-jsonlite::fromJSON(file.path(to,"manifest.json"));stopifnot(mm$spec_hash==fed_hash(comparator_spec()),mm$consensus_id==m$consensus_id)
 z<-readr::read_csv(file.path(to,"coefficients.csv"),show_col_types=FALSE);stopifnot(!any(names(z)%in%c("patient_id","hospitalization_id")))
 rows[[site]]<-z;paths<-c(paths,file.path(to,c("coefficients.csv","manifest.json")))
 }
 pred<-dplyr::bind_rows(rows)
 coef<-pred |> dplyr::group_by(model,term) |> dplyr::summarise(contributing_sites=sum(status=="ok" & is.finite(std_error) & std_error>0),
 num=sum(ifelse(status=="ok" & is.finite(std_error) & std_error>0,estimate/std_error^2,0),na.rm=TRUE),
 den=sum(ifelse(status=="ok" & is.finite(std_error) & std_error>0,1/std_error^2,0),na.rm=TRUE),.groups="drop") |>
 dplyr::mutate(estimate=ifelse(den>0,num/den,0),status=ifelse(contributing_sites==0,"zero_prior_unestimated",ifelse(contributing_sites==1,"single_site","two_or_more_sites"))) |> dplyr::select(model,term,estimate,contributing_sites,status)
 out<-file.path(root,"coordinator",version,"outcomes");fed_write(coef,file.path(out,"comparator_coefficient_support.csv"))
 fed_json(list(predictor_id=substr(fed_hash(list(coef,comparator_spec())),1,20),version=version,training_sites=m$training_sites,held_out_site=m$held_out_site,consensus_id=m$consensus_id,
 coefficients=lapply(split(coef,coef$model),function(z)as.list(setNames(z$estimate,z$term))),input_files=paths,
 comparator_spec=comparator_spec()),file.path(out,"frozen_comparator_predictor.json"))
}
