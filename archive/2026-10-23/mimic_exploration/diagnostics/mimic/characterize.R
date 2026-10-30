source('diagnostics/mimic/portable_lib.R')
root<-mimic_root();a<-readRDS(file.path(root,'local/prepared.rds'));z<-readRDS(file.path(root,'local/derivation.rds'));d<-z$data;r<-a$tables$respiratory_support
catalogue<-r |> dplyr::count(device_category,device_name,name='n');catalogue$n_lower<-10*floor(catalogue$n/10);catalogue$n_upper<-catalogue$n_lower+9;catalogue$n<-NULL
fed_write(catalogue,file.path(root,'aggregate/device_catalogue.csv'))
code<-arrow::open_dataset(file.path(root,'local/source_tables/clif_code_status.parquet')) |> dplyr::select(code_status_name,code_status_category) |> dplyr::distinct() |> dplyr::collect();fed_write(code,file.path(root,'aggregate/code_status_catalogue.csv'))
rows<-list()
for(sub in portable_spec()$subgroups){
 x<-if(sub=='coded_lung_cancer')d else d[d$prior_completed_lung_cancer_encounter,]
 for(v in c('code_at_admission','code_at_landmark','new_dni','new_comfort','code_changed','measured_arf')){
 x$category<-as.character(x[[v]]);x$category[is.na(x$category)]<-'unknown';q<-x |> dplyr::count(phenotype,category,name='n');q$variable<-v;q$subgroup<-sub;q$n_lower<-10*floor(q$n/10);q$n_upper<-q$n_lower+9;q$n<-NULL;rows[[paste(sub,v)]]<-q
 }
}
fed_write(dplyr::bind_rows(rows),file.path(root,'aggregate/code_and_arf_characteristics.csv'))
p<-z$panel;p$location<-'unlocated'
for(id in d$hospitalization_id){
 ix<-which(p$hospitalization_id==id);origin<-d$t0[match(id,d$hospitalization_id)];tt<-origin+(p$h[ix]+.5)*3600
 ad<-a$tables$adt[a$tables$adt$hospitalization_id==id,];ad<-ad[order(ad$in_dttm),]
 for(j in seq_len(nrow(ad))){yes<-which(tt>=ad$in_dttm[j] & (is.na(ad$out_dttm[j])|tt<ad$out_dttm[j]));p$location[ix[yes]]<-tolower(ad$location_category[j])}
}
u<-p |> dplyr::group_by(location) |> dplyr::summarise(hours=dplyr::n(),unknown_hours=sum(support=='UNKNOWN'),unknown_pct=100*mean(support=='UNKNOWN'),.groups='drop')
fed_write(u,file.path(root,'aggregate/unknown_by_location.csv'))
fed_json(list(room_air_source_category_present=any(support_state(r$device_category)=='ROOM_AIR'),room_air_source_name_present=any(grepl('room.?air',r$device_name,ignore.case=TRUE),na.rm=TRUE),room_air_panel_hours=sum(p$support=='ROOM_AIR'),primary_P1_assignments=sum(d$phenotype=='P1'),interpretation='No room-air phenotype coverage is a source/measurement compatibility limitation, not proof that no clinical patients recovered to room air',waterfall_inferred_hours=sum(p$support_carried),no_missing_as_room_air=TRUE),file.path(root,'aggregate/room_air_compatibility.json'))
