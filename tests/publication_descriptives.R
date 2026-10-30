options(clif.aggregate_policy='legacy') # Historical export compatibility test
source('code/lib/load.R')
root<-tempfile('publication_descriptives_');dir.create(root)
make<-function(site,age){
 n<-length(age);t0<-rep(utc('2020-01-02 00:00:00'),n)
 d<-tibble::tibble(patient_id=paste0(site,'p',seq_len(n)),hospitalization_id=paste0(site,'h',seq_len(n)),t0=t0,admission_dttm=t0-3600,discharge_dttm=t0+100*3600,age_years=age,valid_outcome=TRUE,event_type=rep(c('death','live_discharge'),length.out=n),time_from_landmark=rep(c(5,40),length.out=n),phenotype=rep(c('P1','P2'),each=n/2),local_cluster=rep(1:2,each=n/2),sex_category=rep(c('Male','Female'),length.out=n),race_category=rep(c('Black','White'),length.out=n),ethnicity_category='Not Hispanic or Latino',charlson_noncancer=0)
 cfg<-list(site_name=site,analysis=list(run_id='fixture',min_cell=1));p<-list(cohort=d,tables=list(),source_id='fixture')
 panel<-tibble::tibble(hospitalization_id=rep(d$hospitalization_id,each=2),h=rep(0:1,n),support=rep(c('UNKNOWN','IMV'),n),hypox='HU',hyper='CU')
 dest<-file.path(root,site);capture_descriptive_tables(p,d,cfg,dest,panel=panel,phenotypes=c('P1','P2'),consensus_id='shared_fixture')
 target<-file.path(dest,'return/federated/descriptives');dir.create(target,recursive=TRUE);for(f in list.files(file.path(dest,'share/descriptives'),full.names=TRUE))stopifnot(file.copy(f,file.path(target,basename(f))))
 list(data=d,return=file.path(dest,'return'),dest=dest)
}
a<-make('A',20:31);b<-make('B',40:63)
# First-hour support is returned jointly with shared phenotypes, including unknown.
first<-readr::read_csv(file.path(a$dest,'share/descriptives/categorical.csv'),show_col_types=FALSE) |> dplyr::filter(variable=='first_hour_support',cohort=='primary_phenotypes',group!='Overall')
stopifnot(sum(first$n)==nrow(a$data),all(first$category=='UNKNOWN'))
# Phase1 also captures the same cross-tabulation by LOCAL cluster; labels are not pooled.
phase1<-file.path(root,'phase1');p1<-list(cohort=a$data,tables=list(),source_id='fixture')
panel1<-tibble::tibble(hospitalization_id=a$data$hospitalization_id,h=0L,support='HFNC',hypox='HU',hyper='CU')
capture_descriptive_tables(p1,a$data,list(site_name='A',analysis=list(run_id='phase1_fixture',min_cell=1)),phase1,panel=panel1)
z1<-readr::read_csv(file.path(phase1,'share/descriptives/categorical.csv'),show_col_types=FALSE) |> dplyr::filter(variable=='first_hour_support',cohort=='provisional_local_clusters',group!='Overall')
stopifnot(sum(z1$n)==nrow(a$data),all(z1$category=='HFNC'),jsonlite::fromJSON(file.path(phase1,'share/descriptives/manifest.json'))$grouping=='local_cluster')
z<-descriptive_summary(a$data,'primary_phenotypes',c('Overall','P1','P2','I'),'phenotype')
stopifnot(all(z$continuous$observed_n+z$continuous$missing_n==z$continuous$denominator_n))
h<-z$histograms[z$histograms$variable=='age_years'&z$histograms$group=='Overall',];stopifnot(sum(h$n)==nrow(a$data))
# Linked margins must also be masked, not just the rare phenotype cell.
a$data$sofa_total<-1;a$data$sofa_total[1]<-NA
z<-descriptive_protect(descriptive_summary(a$data,'primary_phenotypes',c('Overall','P1','P2'),'phenotype'),3)
stopifnot(all(is.na(z$continuous$mean[z$continuous$variable=='sofa_total'])))
cmd<-c('diagnostics/descriptives/collect.R',shQuote(file.path(root,'combined')),shQuote(c(a$return,b$return)))
rc<-system2(file.path(R.home('bin'),'Rscript'),cmd,stdout=file.path(root,'collect.log'),stderr=file.path(root,'collect.log'));stopifnot(rc==0)
p<-readr::read_csv(file.path(root,'combined/pooled_continuous.csv'),show_col_types=FALSE);p<-p[p$cohort=='landmark_cohort'&p$group=='Overall'&p$variable=='age_years',]
allages<-c(20:31,40:63);stopifnot(abs(p$mean-mean(allages))<1e-10,abs(p$sd-sd(allages))<1e-10,p$observed_n==36)
c<-readr::read_csv(file.path(root,'combined/pooled_categorical.csv'),show_col_types=FALSE);c<-c[c$cohort=='landmark_cohort'&c$group=='Overall'&c$variable=='sex_group',];stopifnot(all(c$denominator_n==36),sum(c$n)==36,sum(c$percent_total)==100)
u<-readr::read_csv(file.path(root,'combined/availability_by_site.csv'),show_col_types=FALSE);stopifnot(all(u$status[u$variable=='stage']=='not_supplied_or_not_derived'))
meta<-jsonlite::fromJSON(file.path(root,'combined/collection_manifest.json'));stopifnot(meta$site_medians_not_averaged,meta$patient_counts_not_assumed_globally_unique)
unlink(root,recursive=TRUE)
cat('PASS: Table1/Table2 moments, denominators, zero groups, linked suppression, missing availability and unequal-site pooled mean/SD\n')
