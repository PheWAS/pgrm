library(pgrm)
library(data.table)


biovu_EUR=annotate_results(results_BioVU_EUR[cohort_match==0], build="hg19",ancestry="EUR")
get_RR(biovu_EUR) # Replicated 651 of 852 for RR=76.4%
get_AER(biovu_EUR) # Expected 1674, replicated 1354 for AE=0.809 (3268 associations for 106 uniq phecodes)
table(biovu_EUR$CI_overlap)
biovu_EUR[is.na(Power)]


biovu_AFR=annotate_results(results_BioVU_AFR[cohort_match==0], build="hg19",ancestry="AFR")
get_RR(biovu_AFR) ## Replicated 11 of 14 for RR=78.6%
get_AER(biovu_AFR) ## Expected 20.2, replicated 19 for AE=0.941 (31 associations for 14 uniq phecodes)
biovu_AFR[is.na(Power)]

anno_BBJ=annotate_results(results_BBJ[cohort_match==0], build="hg19",ancestry="EAS")
get_RR(anno_BBJ) ## Replicated 168 of 219 for RR=76.7%
get_AER(anno_BBJ) ## Expected 261.6, replicated 205 for AE=0.784 (384 associations for 26 uniq phecodes)
anno_BBJ[is.na(Power)]

anno_MGI=annotate_results(results_MGI[cohort_match==0], build="hg38",ancestry="EUR")
get_RR(anno_MGI) ## Replicated 704 of 924 for RR=76.2%
get_AER(anno_MGI) ## Expected 1928.4, replicated 1521 for AE=0.789 (4117 associations for 109 uniq phecodes)
anno_MGI[is.na(Power)]

anno_UKBB=annotate_results(results_UKBB[cohort_match==0], build="hg38",ancestry="EUR")
get_RR(anno_UKBB) # Replicated 706 of 818 for RR=86.3%
get_AER(anno_UKBB) # Expected 1349.5, replicated 1274 for AE=0.944 (2238 associations for 81 uniq phecodes)
anno_UKBB[is.na(Power)]


biovu_EUR$cohort = 'BioVU_EUR'
biovu_AFR$cohort = 'BioVU_AFR'
anno_BBJ$cohort = 'BBJ'
anno_MGI$cohort = 'MGI'
anno_UKBB$cohort = 'UKBB'


d=rbind(biovu_EUR,biovu_AFR,anno_BBJ,anno_MGI,anno_UKBB)
d$cases_needed = NULL
d$cohort_match = NULL
d$SNP = NULL
setcolorder(d, c('assoc_ID', 'rsID', 'risk_allele_dir',
                        'phecode','phecode_string','category_string','ancestry',
                        'cohort','cases','controls','odds_ratio','P','L95','U95',
                 'powered','rep','Power','CI_overlap','rOR','rL95','rU95',
                        'Study_accession','cat_LOG10_P','cat_OR','cat_L95','cat_U95','AF',
                        'pub_count','first_pub_date'))

head(d)
write.table(d,file="annotated_results.csv",sep=",",row.names = F, col.names = T)


## total power assoc
852+14+219+924+844
## total assoc
3268+31+384+4117+2238

rr=c(get_RR(biovu_EUR),get_RR(biovu_AFR),get_RR(anno_BBJ), get_RR(anno_MGI), get_RR(anno_UKBB))
aer=c(get_AER(biovu_EUR),get_AER(biovu_AFR),get_AER(anno_BBJ), get_AER(anno_MGI), get_AER(anno_UKBB))



df <- data.frame(rr, aer)
t.test(rr,aer)
wilcox.test(rr,aer)
cor(rr,aer,method="pearson") ##
cor(rr,aer,method="kendall")
cor(rr,aer,method="spearman")
plot(rr,aer)




