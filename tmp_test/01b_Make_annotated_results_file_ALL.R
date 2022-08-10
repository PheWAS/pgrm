library(pgrm)
library(data.table)

biovu_EUR=annotate_results(results_BioVU_EUR, build="hg19",ancestry="EUR")
biovu_AFR=annotate_results(results_BioVU_AFR, build="hg19",ancestry="AFR")
anno_BBJ=annotate_results(results_BBJ, build="hg19",ancestry="EAS")
anno_MGI=annotate_results(results_MGI, build="hg38",ancestry="EUR")
anno_UKBB=annotate_results(results_UKBB, build="hg38",ancestry="EUR")

biovu_EUR$cohort = 'BioVU_EUR'
biovu_AFR$cohort = 'BioVU_AFR'
anno_BBJ$cohort = 'BBJ'
anno_MGI$cohort = 'MGI'
anno_UKBB$cohort = 'UKBB'
d=rbind(biovu_EUR,biovu_AFR,anno_BBJ,anno_MGI,anno_UKBB)
d$cases_needed = NULL
d$SNP = NULL

setcolorder(d, c('assoc_ID', 'rsID', 'risk_allele_dir',
                 'phecode','phecode_string','category_string','ancestry',
                 'cohort','cases','controls','odds_ratio','P','L95','U95',
                 'powered','rep','Power','CI_overlap','rOR','rL95','rU95',
                 'Study_accession','cat_LOG10_P','cat_OR','cat_L95','cat_U95','RAF',
                 'pub_count','first_pub_date'))

d[,.N,by="cohort"]
d[d$cohort_match==0,.(all=.N,uniq_pheno = length(unique(phecode))),by="cohort"]

head(d)
write.table(d,file="annotated_results_ALL.csv",sep=",",row.names = F, col.names = T)
nrow(d) ## 12497

get_RR(d,include="all") ## Replicated 6219 of 12497 for RR=49.8%

nrow(d[cohort_match==1]) ## 2459
get_RR(d[cohort_match==1],include="all") ## Replicated 1846 of 2459 for RR=75.1%
get_RR(d[cohort_match==0],include="all") ## Replicated 4373 of 10038 for RR=43.6%

table(d[cohort_match==1]$cohort)
prop.table(table(d$cohort,d$cohort_match),1)

nrow(d[cohort_match==0]) ## 10038



