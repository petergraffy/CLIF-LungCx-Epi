source("diagnostics/federated/model_lib.R")
source("code/lib/central.R")
root<-fed_root();version<-commandArgs(trailingOnly=TRUE)[1];out<-file.path(root,"coordinator",version,"outcomes")
manifest<-jsonlite::fromJSON(file.path(root,"coordinator",version,"frozen_manifest.json"));sites<-manifest$training_sites
if(file.exists(file.path(out,"frozen_predictor.json")))stop("Predictor already frozen")
all<-list();prediction<-list();read_log<-character()
for(site in sites){dir<-file.path(root,"exchange","outcomes",version,site)
 m<-jsonlite::fromJSON(file.path(dir,"manifest.json"));stopifnot(m$consensus_id==manifest$consensus_id,m$model_spec_hash==fed_hash(model_spec()))
 all[[site]]<-readr::read_csv(file.path(dir,"associations.csv"),show_col_types=FALSE)
 prediction[[site]]<-readr::read_csv(file.path(dir,"prediction_coefficients.csv"),show_col_types=FALSE)
 read_log<-c(read_log,file.path(dir,c("associations.csv","prediction_coefficients.csv","manifest.json")))
}
a<-dplyr::bind_rows(all);pred<-dplyr::bind_rows(prediction)
pool<-a |> dplyr::group_by(model,endpoint,phenotype,reference,effect_type) |> dplyr::group_modify(~pool_effect(.x)) |> dplyr::ungroup()
coef<-pred |> dplyr::group_by(model,term) |> dplyr::summarise(contributing_sites=sum(status=="ok" & is.finite(std_error) & std_error>0),
 numerator=sum(ifelse(status=="ok" & is.finite(std_error) & std_error>0,estimate/std_error^2,0),na.rm=TRUE),
 denominator=sum(ifelse(status=="ok" & is.finite(std_error) & std_error>0,1/std_error^2,0),na.rm=TRUE),.groups="drop") |>
 dplyr::mutate(estimate=ifelse(denominator>0,numerator/denominator,0),status=ifelse(contributing_sites==0,"zero_prior_unestimated",ifelse(contributing_sites==1,"single_site","two_or_more_sites"))) |>
 dplyr::select(model,term,estimate,contributing_sites,status)
# Zero-prior fallback is used only for unsupported predictive coefficients.
# Association estimates are never fabricated or pooled with <2 contributors.
fed_write(a,file.path(out,"site_associations.csv"));fed_write(pool,file.path(out,"pooled_associations.csv"));fed_write(coef,file.path(out,"prediction_coefficient_support.csv"))
betas<-lapply(split(coef,coef$model),function(z)as.list(stats::setNames(z$estimate,z$term)))
predictor_id<-substr(fed_hash(list(consensus_id=manifest$consensus_id,coefficients=coef,model=model_spec())),1,20)
fed_json(list(schema=FED_SCHEMA,version=version,predictor_id=predictor_id,consensus_id=manifest$consensus_id,training_sites=sites,
 held_out_site=manifest$held_out_site,model_specification=model_spec(),coefficients=betas,input_files=read_log,
 coefficient_support=coef,uncertainty="held-out patient bootstrap conditional on frozen prototypes and predictor"),file.path(out,"frozen_predictor.json"))
message(version,": pooled training-site associations and froze predictor ",predictor_id)
