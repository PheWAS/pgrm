library(pgrm)
library(data.table)


biovu_EUR=annotate_results(results_BioVU_EUR[cohort_match==0], build="hg19",ancestry="EUR")
get_RR(biovu_EUR) # Replicated 651 of 853 for RR=76.3%
get_RR(biovu_EUR,include="all") # Replicated 1354 of 3268 for RR=41.4%
get_AER(biovu_EUR) # Expected 1676.7, replicated 1354 for AE=0.808 (3268 associations for 106 uniq phecodes)
get_powered_rate(biovu_EUR) ## Powered for 853 of 3268 associations 0.26 for 76
table(biovu_EUR$CI_overlap) ##
prop.table(table(biovu_EUR$CI_overlap)) ## % overlap= 0.86811506
biovu_EUR[is.na(Power)]
t.test(biovu_EUR$cat_OR,biovu_EUR$rOR,paired=T) ## mean diff 0.08342065
t.test(biovu_EUR[rep==1]$cat_OR,biovu_EUR[rep==1]$rOR,paired=T) ## mean diff 0.06103702


biovu_AFR=annotate_results(results_BioVU_AFR[cohort_match==0], build="hg19",ancestry="AFR")
get_powered_rate(biovu_AFR) ## 0.483871
get_RR(biovu_AFR) ## Replicated 11 of 14 for RR=78.6%
get_RR(biovu_AFR,include="all") ## Replicated 19 of 31 for RR=61.3%
get_AER(biovu_AFR) ## Expected 20.3, replicated 19 for AE=0.935 (31 associations for 14 uniq phecodes)
biovu_AFR[is.na(Power)]
table(biovu_AFR$CI_overlap) ##  % overlap= 0.7419355
t.test(biovu_AFR$cat_OR,biovu_AFR$rOR,paired=T) ## mean diff 0.8915061
t.test(biovu_AFR[rep==1]$cat_OR,biovu_AFR[rep==1]$rOR,paired=T) ## mean diff 1.245827

anno_BBJ=annotate_results(results_BBJ[cohort_match==0], build="hg19",ancestry="EAS")
get_RR(anno_BBJ) ## Replicated 168 of 219 for RR=76.7%
get_RR(anno_BBJ,include="all") ## Replicated 170 of 222 for RR=76.6%
get_AER(anno_BBJ) ## Expected 262, replicated 205 for AE=0.782 (384 associations for 26 uniq phecodes)
get_powered_rate(anno_BBJ)
table(anno_BBJ$CI_overlap) ##  % overlap= 0.783854167
anno_BBJ[is.na(Power)]
t.test(anno_BBJ$cat_OR,anno_BBJ$rOR,paired=T) ## mean diff 0.2268022
t.test(anno_BBJ[rep==1]$cat_OR,anno_BBJ[rep==1]$rOR,paired=T,conf.level=0.95) ## mean diff 0.08131659

anno_MGI=annotate_results(results_MGI[cohort_match==0], build="hg38",ancestry="EUR")
get_RR(anno_MGI) ## Replicated 706 of 928 for RR=76.1%
get_RR(anno_MGI,include="all") ## Replicated 1521 of 4117 for RR=36.9%
get_AER(anno_MGI) ## Expected 1932.3, replicated 1521 for AE=0.787 (4117 associations for 109 uniq phecodes)
get_powered_rate(anno_MGI)
table(anno_MGI$CI_overlap)
anno_MGI[is.na(Power)]
t.test(anno_MGI$cat_OR,anno_MGI$rOR,paired=T) ## mean diff 0.07281874
t.test(anno_MGI[rep==1]$cat_OR,anno_MGI[rep==1]$rOR,paired=T) ## mean diff 0.03712937
prop.table(table(anno_MGI$CI_overlap)) ##  % overlap= 0.882681564


anno_UKBB=annotate_results(results_UKBB[cohort_match==0], build="hg38",ancestry="EUR")
get_RR(anno_UKBB) # Replicated 727 of 853 for RR=85.2%
get_RR(anno_UKBB,include="all") ## Replicated 1274 of 2238 for RR=56.9%
get_AER(anno_UKBB) #Expected 1353.3, replicated 1274 for AE=0.941 (2238 associations for 81 uniq phecodes)
get_powered_rate(anno_UKBB)
table(anno_UKBB$CI_overlap)
t.test(anno_UKBB$cat_OR,anno_UKBB$rOR,paired=T) ## mean diff 0.06598557
t.test(anno_UKBB[rep==1]$cat_OR,anno_UKBB[rep==1]$rOR,paired=T) ## mean diff 0.03162266
anno_UKBB[is.na(Power)]
head(anno_UKBB)

t.test(anno_UKBB$cat_OR,anno_UKBB$rOR,paired=T) ## mean diff 0.06598557
t.test(anno_UKBB[rep==1]$cat_OR,anno_UKBB[rep==1]$rOR,paired=T) ## mean diff 0.03162266


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
                        'Study_accession','cat_LOG10_P','cat_OR','cat_L95','cat_U95','RAF',
                        'pub_count','first_pub_date'))

head(d)
write.table(d,file="annotated_results.csv",sep=",",row.names = F, col.names = T)


## total power assoc
852+14+219+924+844
## total assoc
3268+31+384+4117+2238

rr=c(get_RR(biovu_EUR),get_RR(biovu_AFR),get_RR(anno_BBJ), get_RR(anno_MGI), get_RR(anno_UKBB))
aer=c(get_AER(biovu_EUR),get_AER(biovu_AFR),get_AER(anno_BBJ), get_AER(anno_MGI), get_AER(anno_UKBB))


head(d)


df <- data.frame(rr, aer)
t.test(rr,aer,paired=T)
wilcox.test(rr,aer,paired=T)
cor(rr,aer,method="pearson") ##
cor(rr,aer,method="kendall")
cor(rr,aer,method="spearman")
plot(rr,aer)

plot(d$cat_OR,d$odds_ratio)

d$cat_beta=log(d$cat_OR)
d[risk_allele_dir=="ref"]$cat_beta= d[risk_allele_dir=="ref"]$cat_beta*-1
d$test_beta = log(d$odds_ratio)
#d[risk_allele_dir=="ref"]$test_beta= d[risk_allele_dir=="ref"]$test_beta*-1
hist(d$test_beta)
hist(d$cat_beta)

plot(d[P<0.05]$cat_OR,d[P<0.05]$test_OR)
Model=lm(d[P<0.05]$cat_OR~d[P<0.05]$test_OR)
abline(Model)
legend("topleft",legend=paste("R2 is", format(summary(Model)$r.squared,digits=3)))

plot(d$cat_OR,d$test_OR)
Model=lm(d$cat_OR~d$test_OR)
abline(Model)
legend("topleft",legend=paste("R2 is", format(summary(Model)$r.squared,digits=3)))

t=t.test(d[P<0.05]$cat_beta,d[P<0.05]$test_beta,paired=T)
t

t=t.test(d[P<0.05]$cat_OR,d[P<0.05]$test_OR,paired=T)
t



t=t.test(d$cat_OR,d$test_OR,paired=T)
t
t$p.value

cor(d$cat_beta,d$test_beta)
d$test_OR = d$odds_ratio
d[risk_allele_dir=="ref"]$test_OR = 1/d[risk_allele_dir=="ref"]$test_OR
t=t.test(d$cat_OR,d$test_OR,paired=T)
t$p.value
t
cor(d$cat_OR,d$test_OR,method="kendall") ##

m=glm(test_OR~cat_OR,data=d)
summary(m)
library(effects)
plot(allEffects(m))




