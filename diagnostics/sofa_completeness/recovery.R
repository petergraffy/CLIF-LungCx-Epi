# Conservative measurement recovery; no missing-domain-zero or oxygen assumptions.
# Exact-timestamp BP/GCS tuples must contain one distinct valid value per field.
sofa_derived_rows <- function(x, category, value, components, output, ranges, calculate) {
 z <- x[x[[category]] %in% components, ,drop=FALSE]
 z$value_recovery <- num(z[[value]])
 for(j in seq_along(components)){
  bad <- z[[category]]==components[j] & (z$value_recovery<ranges[[j]][1] | z$value_recovery>ranges[[j]][2])
  z$value_recovery[which(bad)] <- NA_real_
 }
 z <- z[is.finite(z$value_recovery),,drop=FALSE]
 z <- z |> dplyr::group_by(hospitalization_id,recorded_dttm,.data[[category]]) |>
  dplyr::summarise(n=dplyr::n_distinct(value_recovery),v=dplyr::first(value_recovery),.groups='drop') |>
  dplyr::filter(n==1) |> dplyr::select(-n) |>
  tidyr::pivot_wider(names_from=dplyr::all_of(category),values_from=v)
 for(c in components)if(!c%in%names(z))z[[c]]<-NA_real_
 keep <- stats::complete.cases(z[,components,drop=FALSE]);z<-z[keep,,drop=FALSE]
 z[[value]] <- calculate(z);z[[category]]<-output
 z <- z[is.finite(z[[value]]),c('hospitalization_id','recorded_dttm',category,value),drop=FALSE]
 tibble::as_tibble(z)
}
sofa_recover_tuples <- function(tables) {
 v<-tables$vitals;g<-tables$patient_assessments
 mp<-sofa_derived_rows(v,'vital_category','vital_value',c('sbp','dbp'),'map',list(c(40,300),c(20,200)),function(z)ifelse(z$sbp>z$dbp,(z$sbp+2*z$dbp)/3,NA_real_))
 gc<-sofa_derived_rows(g,'assessment_category','numerical_value',c('gcs_eye','gcs_verbal','gcs_motor'),'gcs_total',list(c(1,4),c(1,5),c(1,6)),function(z)ifelse(z$gcs_eye%%1==0 & z$gcs_verbal%%1==0 & z$gcs_motor%%1==0,z$gcs_eye+z$gcs_verbal+z$gcs_motor,NA_real_))
 # Only fill a domain when the original window contains no valid observation.
 list(map=mp,gcs=gc)
}
