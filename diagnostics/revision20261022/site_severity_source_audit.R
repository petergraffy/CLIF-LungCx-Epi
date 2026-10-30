source('diagnostics/revision20261022/lib.R')
site<-commandArgs(trailingOnly=TRUE)[1];home<-file.path(rev_root(),'sites',site)
a<-readRDS(file.path(old_root(),'sites',site,'local/discovery.rds'));d<-a$data;m<-jsonlite::fromJSON(file.path(a$input_root,'local/input_manifest.json'))
source_dir<-Sys.getenv(paste0('CLIF_SITE_TABLES_',site),unset=dirname(m$source_file))
x<-arrow::open_dataset(file.path(source_dir,'clif_labs.parquet'))
labels<-x |> dplyr::select(lab_name,lab_category,lab_loinc_code) |> dplyr::distinct() |> dplyr::collect()
labels<-labels[grepl('hct|hematocrit|haematocrit',labels$lab_name,ignore.case=TRUE)|grepl('4544-3|20570-8',labels$lab_loinc_code),]
n<-0
if(nrow(labels)){
 names<-unique(labels$lab_name);ids<-d$hospitalization_id
 z<-x |> dplyr::filter(lab_name%in%names,hospitalization_id%in%ids) |> dplyr::select(hospitalization_id,lab_result_dttm,lab_value_numeric) |> dplyr::collect() |> dplyr::inner_join(d |> dplyr::select(hospitalization_id,t0),by='hospitalization_id')
 z$lab_result_dttm<-utc(z$lab_result_dttm)
 n<-length(unique(z$hospitalization_id[which(z$lab_result_dttm>=z$t0 & z$lab_result_dttm<z$t0+24*3600 & is.finite(num(z$lab_value_numeric)))]))
}
b<-fed_band_n(n)
fed_write(tibble::tibble(site=site,component='hematocrit_source_label_or_LOINC_audit',available_lower=b[1],available_upper=b[2],status=if(n==0)'no_first24h_source_measurements_found' else if(site=='NU')'arterial_hematocrit_present_under_other_category' else 'canonical_hematocrit_present'),file.path(home,'export/severity_source_audit.csv'))
