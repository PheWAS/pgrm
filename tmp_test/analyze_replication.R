library(ggplot2)

library(pgrm)
results_MGI
results_UKBB
results_BioVU_EUR

PGRM=get_PGRM(ancestry='EUR',build='hg19')
MGI=annotate_results(results_MGI[cohort_match==0], build='hg38',ancestry='EUR')
UKBB=annotate_results(results_UKBB[cohort_match==0], build='hg38',ancestry='EUR')
BioVU_EUR=annotate_results(results_BioVU_EUR[cohort_match==0], build='hg19',ancestry='EUR')

get_RR(UKBB)

foo=annotate_results(results_UKBB[cohort_match==0 & SNP %in% c('11:113936459:G:A','2:8311368:A:G')], build='hg38',ancestry='EUR')

foo

BioVU_AFR=annotate_results(results_BioVU_AFR[cohort_match==0], build='hg19',ancestry='AFR')
BBJ=annotate_results(BBJ[cohort_match==0], build='hg19',ancestry='EAS')



cols_to_keep=
  c('assoc_ID','SNP','phecode','phecode_string','category_string','dataset','P','odds_ratio','L95','U95')

PGRM$P = 10^PGRM$cat_LOG10_P
names(PGRM)[12]='odds_ratio'
names(PGRM)[13]='L95'
names(PGRM)[14]='U95'
PGRM$dataset='catalog'

MGI$dataset='MGI'
UKBB$dataset='UKBB'
BioVU$dataset='BioVU'

foo=c("SNP","phecode")

d=rbind(PGRM[,c('assoc_ID','SNP','phecode','phecode_string','category_string','dataset','P','odds_ratio','L95','U95')],
        MGI[,c('assoc_ID','SNP','phecode','phecode_string','category_string','dataset','P','odds_ratio','L95','U95')])
d=rbind(d,UKBB[,c('assoc_ID','SNP','phecode','phecode_string','category_string','dataset','P','odds_ratio','L95','U95')])
d=rbind(d,BioVU[,c('assoc_ID','SNP','phecode','phecode_string','category_string','dataset','P','odds_ratio','L95','U95')])
table(d$dataset,d$phecode)
head(d)
