library(pgrm)
library(data.table)


biovu_EUR=annotate_results(results_BioVU_EUR[cohort_match==0], build="hg19",ancestry="EUR")
get_RR(biovu_EUR) # Replicated 657 of 865 for RR=76.0%
get_AER(biovu_EUR) # Expected 1689.8, replicated 1354 for AE=0.801 (3268 associations for 106 uniq phecodes)
table(biovu_EUR$CI_overlap)
biovu_EUR[is.na(Power)]


biovu_AFR=annotate_results(results_BioVU_AFR[cohort_match==0], build="hg19",ancestry="AFR")
get_RR(biovu_AFR) ## Replicated 11 of 15 for RR=73.3%
get_AER(biovu_AFR) ## Expected 20.3, replicated 19 for AE=0.935 (31 associations for 14 uniq phecodes)
biovu_AFR[is.na(Power)]

anno_BBJ=annotate_results(results_BBJ[cohort_match==0], build="hg19",ancestry="EAS")
get_RR(anno_BBJ) ## Replicated 170 of 221 for RR=76.9%
get_AER(anno_BBJ) ## Expected 264, replicated 205 for AE=0.777 (384 associations for 26 uniq phecodes)
anno_BBJ[is.na(Power)]

anno_MGI=annotate_results(results_MGI[cohort_match==0], build="hg38",ancestry="EUR")
get_RR(anno_MGI) ## Replicated 708 of 938 for RR=75.5%
get_AER(anno_MGI) ## Expected 1947.7, replicated 1521 for AE=0.781 (4117 associations for 109 uniq phecodes)
anno_MGI[is.na(Power)]

anno_UKBB=annotate_results(results_UKBB[cohort_match==0], build="hg38",ancestry="EUR")
get_RR(anno_UKBB) # Replicated 739 of 874 for RR=84.6%
get_AER(anno_UKBB) # Expected 1368, replicated 1274 for AE=0.931 (2238 associations for 81 uniq phecodes)
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







annotate_power(anno_MGI[is.na(Power)],max_thresh=50)

annotate_power = function(annotated_results, LOUD = FALSE,max_thresh=100) {
  Power = NULL

  annotated_results[, Power := as.numeric(NA)]
  total = nrow(annotated_results)
  if (isTRUE(LOUD)) {
    print('Doing power calculations')}
  for (i in seq_len(total)) {
    if (isTRUE(LOUD) && i %% 100 == 0) {
      print(glue('{i} of {total}'))}

    annotated_results_tmp = annotated_results[i]
    odds_ratio = annotated_results_tmp$cat_L95
    AF = annotated_results_tmp$RAF
    k = annotated_results_tmp$controls / annotated_results_tmp$cases

    ## flip AF and OR to minor allele
    if (AF > 0.5) {
      AF = 1 - AF
      odds_ratio = 1 / odds_ratio
    }

    N = annotated_results_tmp$controls + annotated_results_tmp$cases
    ## control:case ratio ceiling of max_thresh
    if (k > max_thresh) {
      k = max_thresh
      N = annotated_results_tmp$cases * max_thresh
    }

    pwr = genpwr.calc(calc = 'power', model = 'logistic', ge.interaction = NULL,
                      Case.Rate = NULL, k = k, N = N, MAF = AF, OR = odds_ratio,
                      Alpha = 0.05, Power = NULL, True.Model = c('Additive'),
                      Test.Model = c('Additive'))
    annotated_results[i]$Power = pwr$Power_at_Alpha_0.05}
  return(annotated_results)}


annotate_power = function(annotated_results, LOUD = FALSE,max_thresh=100) {
  Power = NULL

  annotated_results[, Power := as.numeric(NA)]
  total = nrow(annotated_results)
  if (isTRUE(LOUD)) {
    print('Doing power calculations')}
  for (i in seq_len(total)) {
    if (isTRUE(LOUD) && i %% 100 == 0) {
      print(glue('{i} of {total}'))}

    annotated_results_tmp = annotated_results[i]
    odds_ratio = annotated_results_tmp$cat_L95
    AF = annotated_results_tmp$RAF
    ## flip AF and OR to minor allele
    if (AF > 0.5) {
      AF = 1 - AF
      odds_ratio = 1 / odds_ratio
    }
    k = annotated_results_tmp$controls / annotated_results_tmp$cases
    N = annotated_results_tmp$controls + annotated_results_tmp$cases
    ## control:case ratio ceiling of max_thresh
    if (k > max_thresh) {
      k = max_thresh
      N = annotated_results_tmp$cases * max_thresh
    }

    pwr = genpwr.calc(calc = 'power', model = 'logistic', ge.interaction = NULL,
                      Case.Rate = NULL, k = k, N = N, MAF = AF, OR = odds_ratio,
                      Alpha = 0.05, Power = NULL, True.Model = c('Additive'),
                      Test.Model = c('Additive'))
    annotated_results[i]$Power = pwr$Power_at_Alpha_0.05
    #  assoc_ID = annotated_results_tmp$assoc_ID
    #  x=pwr$Power_at_Alpha_0.05
    #  print(glue('{assoc_ID} {x} k={k} N={N} AF={AF} OR={odds_ratio}'))

  }
  return(annotated_results)}

print(glue('Replicated {numerator} of {denominator} for RR={f_RR}',
           f_RR = sprintf('%1.1f%%', 100 * RR)))}


results1=biovu_EUR
results2=anno_UKBB
compare_annotated_results(biovu_EUR, anno_UKBB)
compare_annotated_results = function(results1, results2){

  summary=data.table()
  results1$dataset='results1'
  results2$dataset='results2'
  r_long=rbind(results1[,c('assoc_ID','odds_ratio','P','L95','U95','rep','powered','Power','rOR','rL95','rU95','dataset')],
               results2[,c('assoc_ID','odds_ratio','P','L95','U95','rep','powered','Power','rOR','rL95','rU95','dataset')])
  r_long=data.table(r_long)
 ## Compare RR all
  m=glm(rep~dataset,data=r_long,family="binomial")
  P=summary(m)$coeff['datasetresults2','Pr(>|z|)']
  OR=exp(summary(m)$coeff['datasetresults2','Estimate'])
  CIs=exp(confint.default(m))
  L95=CIs['datasetresults2',1]
  U95=CIs['datasetresults2',2]
  print(glue('\n--------------------------------------'))
  print(glue('Replication rate (all) comparison'))
  RR_ALL_r1=get_RR(results1,include="all",LOUD=FALSE)
  RR_ALL_r2=get_RR(results2,include="all",LOUD=FALSE)
  print(glue('Result1 replication rate (overall) = {r1_RR}',
             r1_RR = sprintf('%1.1f%%', 100 * RR_ALL_r1)))
  print(glue('Result2 replication rate (overall) = {r2_RR}',
             r2_RR = sprintf('%1.1f%%', 100 * RR_ALL_r2)))
  print(glue('Logistic regression rep~dataset'))
  print(glue('Odds ratio (95% CIs) {round(OR,4)} ({round(L95,4)} to {round(U95,4)})'))
  print(glue('P-value {P}'))

  ## Compare Power
  m=glm(powered~dataset,data=r_long,family="binomial")
  P=summary(m)$coeff['datasetresults2','Pr(>|z|)']
  OR=exp(summary(m)$coeff['datasetresults2','Estimate'])
  CIs=exp(confint.default(m))
  L95=CIs['datasetresults2',1]
  U95=CIs['datasetresults2',2]
  print(glue('\n--------------------------------------'))
  print(glue('Powered comparison'))
  Power_r1=get_powered_rate(results1,LOUD=FALSE)
  Power_r2=get_powered_rate(results2,LOUD=FALSE)
  print(glue('Result1 replication rate (overall) = {r1_power}',
             r1_power = sprintf('%1.1f%%', 100 * Power_r1)))
  print(glue('Result2 replication rate (overall) = {r2_power}',
             r2_power = sprintf('%1.1f%%', 100 * Power_r2)))
  print(glue('Odds ratio (95% CIs) {round(OR,4)} ({round(L95,4)} to {round(U95,4)})'))
  print(glue('P-value {P}'))

  ## Compare RR powered
  m=glm(rep~dataset,data=r_long[powered==1],family="binomial")
  P=summary(m)$coeff['datasetresults2','Pr(>|z|)']
  OR=exp(summary(m)$coeff['datasetresults2','Estimate'])
  CIs=exp(confint.default(m))
  L95=CIs['datasetresults2',1]
  U95=CIs['datasetresults2',2]
  print(glue('\n--------------------------------------'))
  print(glue('Replication rate (all) comparison'))
  RR_ALL_r1=get_RR(results1,LOUD=FALSE)
  RR_ALL_r2=get_RR(results2,LOUD=FALSE)
  print(glue('Result1 replication rate (powered) = {r1_RR}',
             r1_RR = sprintf('%1.1f%%', 100 * RR_ALL_r1)))
  print(glue('Result2 replication rate (powered) = {r2_RR}',
             r2_RR = sprintf('%1.1f%%', 100 * RR_ALL_r2)))
  print(glue('Odds ratio (95% CIs) {round(OR,4)} ({round(L95,4)} to {round(U95,4)})'))
  print(glue('P-value {P}'))

  ## Compare RR, control for power
  m=glm(rep~dataset+Power,data=r_long,family="binomial")
  P=summary(m)$coeff['datasetresults2','Pr(>|z|)']
  OR=exp(summary(m)$coeff['datasetresults2','Estimate'])
  CIs=exp(confint.default(m))
  L95=CIs['datasetresults2',1]
  U95=CIs['datasetresults2',2]
  print(glue('\n--------------------------------------'))
  print(glue('Replication comparison, controlling for Power'))
  AER_r1=get_AER(results1,LOUD=FALSE)
  AER_r2=get_AER(results2,LOUD=FALSE)
  print(glue('Result1 Actual:Expected = {round(r1_AER,3)}'))
  print(glue('Result2 Actual:Expected = {round(r2_AER,3)}'))
  print(glue('Odds ratio (95% CIs) {round(OR,4)} ({round(L95,4)} to {round(U95,4)})'))
  print(glue('P-value {P}'))

  r1=results1[,c('assoc_ID','rsID','phecode','phecode_string','category_string','odds_ratio','P','L95','U95','rep','powered','Power','rOR','rL95','rU95')]
  r2=results2[,c('assoc_ID','odds_ratio','P','L95','U95','rep','powered','Power','rOR','rL95','rU95')]
  r=merge(r1,r2,by="assoc_ID")
  overlapping_rows=nrow(r)
  if(overlapping_rows==0){
    print("no overlapping rows in datasets")
    return(0)
  }
  both_powered=nrow(r[powered.x==1 & powered.y==1])
  both_rep=nrow(r[rep.x==1 & rep.y==1])
  both_powered_rep=nrow(r[rep.x==1 & rep.y==1 & powered.x==1 & powered.y==1])

  r$both_powered = 0
  r[powered.x==1 & powered.y==1]$both_powered = 1

  r$both_rep = 0
  r[rep.x==1 & rep.y==1]$both_rep = 1

  r$both_powered_rep = 0
  r[rep.x==1 & rep.y==1 & powered.x==1 & powered.y==1]$both_powered_rep = 1

  table(r$both_powered,r$both_rep)

  print(glue('\n--------------------------------------'))
  print(glue('Compare odds ratios, all'))

  t=t.test(r$rOR.x,r$rOR.y,paired=T)
  print(glue('Odds ratio comparison (all), n={overlapping_rows}'))
  print(glue('P-value {t$p.value}'))
  print(glue('Mean difference {round(t$estimate,4)}'))
  print(glue('95% CI {round(t$conf.int[1],4)} to {round(t$conf.int[2],4)}'))

  print(glue('\n--------------------------------------'))
  print(glue('Compare odds ratios, Powered'))

  t=t.test(r[rep.x==1 & rep.y==1]$rOR.x,r[rep.x==1 & rep.y==1]$rOR.y,paired=T)
  print(glue('Odds ratio comparison (both replicated), n={both_rep}'))
  print(glue('P-value {t$p.value}'))
  print(glue('Mean difference {round(t$estimate,4)}'))
  print(glue('95% CI {round(t$conf.int[1],4)} to {round(t$conf.int[2],4)}'))}
