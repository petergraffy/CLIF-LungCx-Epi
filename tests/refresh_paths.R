source('tests/revision_gaps.R')
local({
 vars<-c('CLIF_REVISION_PRIMARY_INPUT','CLIF_FEDERATED_INPUT_UCMC');saved<-Sys.getenv(vars,unset=NA_character_)
 on.exit({for(i in seq_along(vars))if(is.na(saved[i]))Sys.unsetenv(vars[i]) else do.call(Sys.setenv,setNames(list(saved[i]),vars[i]))},add=TRUE)
 Sys.setenv(CLIF_REVISION_PRIMARY_INPUT=tempdir(),CLIF_FEDERATED_INPUT_UCMC=tempdir())
 stopifnot(old_root()==normalizePath(tempdir()),revision_site_input('UCMC')==normalizePath(tempdir()))
})
d$phenotype[1]<-'P5'
x<-score_matrix(d,'summaries_phenotype','missing_domain_zero_sensitivity',paste0('P',1:5))
stopifnot('phenotype_P5'%in%colnames(x),x[1,'phenotype_P5']==1)
x3<-score_matrix(d,'summaries_phenotype','missing_domain_zero_sensitivity',paste0('P',1:3))
stopifnot(!any(c('phenotype_P4','phenotype_P5')%in%colnames(x3)))
cat('Refresh source-path overrides and variable phenotype-count contracts passed\n')
