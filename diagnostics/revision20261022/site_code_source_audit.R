source('diagnostics/revision20261022/lib.R')
site<-commandArgs(trailingOnly=TRUE)[1];home<-file.path(rev_root(),'sites',site);m<-jsonlite::fromJSON(file.path(revision_site_input(site),'local/input_manifest.json'))
source_dir<-Sys.getenv(paste0('CLIF_SITE_TABLES_',site),unset=dirname(m$source_file))
files<-list.files(source_dir,pattern='[.]parquet$',recursive=TRUE,full.names=TRUE,ignore.case=TRUE)
bases<-tolower(sub('^clif_','',tools::file_path_sans_ext(basename(files)),ignore.case=TRUE));files<-files[bases=='code_status']
if(length(files)>1)stop('Expected at most one code_status table.')
source_available<-length(files)==1
if(source_available){
 x<-arrow::open_dataset(files);available<-names(x)
 z<-x |> dplyr::select(dplyr::any_of(c('code_status_name','code_status_category'))) |> dplyr::distinct() |> dplyr::collect()
 z<-with_columns(z,c('code_status_name','code_status_category'))
}else{available<-character();z<-tibble::tibble(code_status_name=character(),code_status_category=character())}
name_available<-'code_status_name'%in%available
z$site<-rep(site,nrow(z));raw<-toupper(gsub('[^A-Za-z]','',z$code_status_category));nm<-toupper(z$code_status_name)
z$interpretation<-if(!name_available)rep('source_order_names_unavailable',nrow(z)) else ifelse(raw%in%c('DNRDNI','DNARDNI','DNI') & !grepl('DNI|DO NOT INTUBATE',nm),'canonical_DNI_not_explicit_in_source_order_name',ifelse(raw=='OTHER','limitations_not_resolved_by_canonical_other','no_label_conflict_identified'))
fed_write(z,file.path(home,'export/code_source_catalogue.csv'))
fed_json(list(site=site,source_available=source_available,source_order_names_available=name_available,comfort_category_present=if(source_available)any(raw%in%c('AND','ALLOWNATURALDEATH','CMO','COMFORTMEASURESONLY')) else NA,explicit_DNI_in_order_name=if(name_available)any(grepl('DNI|DO NOT INTUBATE',nm),na.rm=TRUE) else NA,canonical_DNI_with_nonspecific_order_name=if(name_available)any(z$interpretation=='canonical_DNI_not_explicit_in_source_order_name',na.rm=TRUE) else NA,interpretation='Catalogue audit, not source EHR validation; unavailable source is not absence of a preference.'),file.path(home,'export/code_source_audit.json'))
