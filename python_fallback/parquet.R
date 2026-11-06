clif_fallback_python<-function(){p<-Sys.getenv('CLIF_FALLBACK_PYTHON');if(!nzchar(p)||!file.exists(p))stop('Use the python_fallback/site.R entry point with a configured Python executable.');p}
# PyArrow Parquet transport. R Arrow is neither loaded nor required.
clif_python_bridge <- function(request, python) {
 dir<-tempfile('clif_pyarrow_');dir.create(dir,mode='0700')
 on.exit(unlink(dir,recursive=TRUE))
 request$metadata<-file.path(dir,'schema.json');request$csv<-file.path(dir,'rows.csv')
 request$numeric_dir<-dir
 request$null<-paste0('__CLIF_NULL_',basename(dir),'__')
 if(!is.null(request$ids)){
  request$ids_file<-file.path(dir,'ids.json')
  jsonlite::write_json(as.character(request$ids),request$ids_file,auto_unbox=FALSE,na='null')
  request$ids<-NULL
 }
 if(!is.null(request$data)){
  x<-request$data;request$data<-NULL
  fields<-lapply(names(x),function(n){v<-x[[n]]
   kind<-if(inherits(v,'POSIXt'))'timestamp' else if(inherits(v,'Date'))'date' else if(is.logical(v))'boolean' else if(is.integer(v)||inherits(v,'integer64'))'integer' else if(is.numeric(v))'double' else 'string'
   list(name=n,kind=kind,timezone=if(kind=='timestamp')'UTC' else NULL)
  })
  jsonlite::write_json(fields,request$metadata,auto_unbox=TRUE,null='null')
  for(n in names(x)){
   if(inherits(x[[n]],'POSIXt'))x[[n]]<-as.numeric(x[[n]])
   if(inherits(x[[n]],'integer64'))x[[n]]<-as.character(x[[n]])
   if(is.double(x[[n]])&&!inherits(x[[n]],'Date')){
    v<-x[[n]];i<-match(n,names(x))-1L
    writeBin(v,file.path(dir,paste0(i,'.bin')),size=8,endian='little')
    writeBin(as.raw(is.na(v)&!is.nan(v)),file.path(dir,paste0(i,'.null')))
    x[[n]]<-rep(NA_character_,length(v))
   }
   if(is.character(x[[n]])&&any(x[[n]]==request$null,na.rm=TRUE))stop('Null marker collision')
  }
  # Numeric values travel as binary doubles; CSV carries text/integer fields.
  old<-options(digits=17,scipen=999);on.exit(options(old),add=TRUE)
  utils::write.table(x,request$csv,sep=',',row.names=FALSE,col.names=TRUE,na=request$null,qmethod='double',fileEncoding='UTF-8')
 }
 jsonlite::write_json(request,file.path(dir,'request.json'),auto_unbox=TRUE,null='null',digits=17)
 status<-system2(python,c(shQuote('fallback/parquet_bridge.py'),shQuote(file.path(dir,'request.json'))),stdout=file.path(dir,'worker.log'),stderr=file.path(dir,'worker.log'))
 if(status!=0)stop('PyArrow bridge failed: ',paste(readLines(file.path(dir,'worker.log')),collapse='\n'))
 if(request$operation=='write')return(invisible(request$path))
 fields<-jsonlite::fromJSON(request$metadata,simplifyVector=FALSE)
 if(isTRUE(request$schema_only))return(fields)
 parsers<-setNames(lapply(fields,function(f){
  switch(f$kind,timestamp=readr::col_character(),date=readr::col_date(),boolean=readr::col_logical(),double=readr::col_character(),integer=if(f$bits<=32)readr::col_integer() else readr::col_character(),readr::col_character())
 }),vapply(fields,`[[`,character(1),'name'))
 x<-as.data.frame(readr::read_csv(request$csv,col_types=do.call(readr::cols,parsers),na=request$null,show_col_types=FALSE,progress=FALSE,trim_ws=FALSE,num_threads=1))
 if(nrow(readr::problems(x)))stop('CSV transport parsing problems')
 for(f in fields){n<-f$name
  if(f$kind%in%c('timestamp','double')){
   i<-match(n,names(x))-1L
   x[[n]]<-readBin(file.path(dir,paste0(i,'.bin')),'double',n=nrow(x),size=8,endian='little')
   missing<-readBin(file.path(dir,paste0(i,'.null')),'raw',n=nrow(x))
   stopifnot(length(x[[n]])==nrow(x),length(missing)==nrow(x))
   x[[n]][as.logical(missing)]<-NA_real_
  }
  if(f$kind=='timestamp'){
   x[[n]]<-as.POSIXct(x[[n]],origin='1970-01-01',tz=if(is.null(f$timezone))'UTC' else f$timezone)
   # Match R Arrow's representation of timezone-naive columns. The existing
   # source_time() policy must receive the same schema AND POSIX representation.
   if(is.null(f$timezone))attr(x[[n]],'tzone')<-NULL
  }
  if(f$kind=='integer'&&f$bits>32){
   # Match Arrow's default: represent int64 as doubles only when all are exact.
   z<-bit64::as.integer64(x[[n]])
   if(any(!is.na(x[[n]])&is.na(z)))stop('Integer value outside the supported signed int64 range: ',n)
   if(all(is.na(z)|(z<=bit64::as.integer64('9007199254740992')&z>=bit64::as.integer64('-9007199254740992'))))x[[n]]<-as.numeric(z) else x[[n]]<-z
  }
 }
 attr(x,'pyarrow_fields')<-fields
 x
}
clif_parquet_schema <- function(path,python=clif_fallback_python()){
 fields<-clif_python_bridge(list(operation='read',path=path,schema_only=TRUE),python)
 make_type<-function(s)local({value<-s;list(ToString=function()value)})
 list(names=vapply(fields,`[[`,character(1),'name'),fields=fields,
  GetFieldByName=function(n)list(type=make_type(fields[[match(n,vapply(fields,`[[`,character(1),'name'))]]$type)))
}
clif_read_parquet <- function(path,python=clif_fallback_python(),columns=NULL,id_column=NULL,ids=NULL,distinct=FALSE){
 x<-clif_python_bridge(list(operation='read',path=path,columns=if(is.null(columns))NULL else as.list(columns),id_column=id_column,ids=ids,distinct=distinct),python);attr(x,'pyarrow_fields')<-NULL;tibble::as_tibble(x)
}
clif_write_parquet <- function(x,path,python=clif_fallback_python())clif_python_bridge(list(operation='write',path=path,data=x),python)

