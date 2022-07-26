
library(ggplot2)

library(pgrm)
results_MGI
results_UKBB
results_BioVU_EUR

PGRM=get_PGRM(ancestry='EUR',build='hg19')
MGI=annotate_results(results_MGI[cohort_match==0], build='hg38',ancestry='EUR')
UKBB=annotate_results(results_UKBB[cohort_match==0], build='hg38',ancestry='EUR')
BioVU=annotate_results(results_BioVU_EUR[cohort_match==0], build='hg19',ancestry='EUR')



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

cur_phecode='010'

ggplot(d[phecode==cur_phecode], aes(x=odds_ratio,y=SNP,label=P,color=dataset))+geom_point()

ggplot(d[phecode==cur_phecode], aes(x=odds_ratio,y=SNP,label=P,color=dataset))+
  geom_errorbar(aes(xmin=L95), xmax=U95), width=.1)

  geom_point()+geom_text(hjust=-.4, vjust=-1)+ylab('')+xlab('Odds ratio')+
  scale_y_discrete(limits=rev)+theme_classic()+
  geom_vline(xintercept=1,color="grey",linetype="dashed")



head(PGRM)
head(MGI)
names(MGI)[3]='MGI_cases'
names(MGI)[4]='MGI_controls'
names(MGI)[5]='MGI_odds_ratio'
names(MGI)[6]='MGI_P'
names(MGI)[7]='MGI_L95'
names(MGI)[8]='MGI_U95'
names(MGI)[26]='MGI_Power'
head(MGI)
d=merge(PGRM[ancestry=='EUR',c('assoc_ID','SNP','phecode','phecode_string',
                               'category_string','cat_LOG10_P','cat_OR','cat_L95','cat_U95')],
        MGI[,c('assoc_ID','MGI_cases','MGI_controls','MGI_odds_ratio','MGI_P','MGI_L95','MGI_Power')], by='assoc_ID', all.x=T)

names(UKBB)[3]='UKBB_cases'
names(UKBB)[4]='UKBB_controls'
names(UKBB)[5]='UKBB_odds_ratio'
names(UKBB)[6]='UKBB_P'
names(UKBB)[7]='UKBB_L95'
names(UKBB)[8]='UKBB_U95'
names(UKBB)[26]='UKBB_Power'

d=merge(d,UKBB[,c('assoc_ID','UKBB_cases','UKBB_controls','UKBB_odds_ratio','UKBB_P','UKBB_L95','UKBB_Power')], by='assoc_ID', all.x=T)

names(BioVU)[3]='BioVU_cases'
names(BioVU)[4]='BioVU_controls'
names(BioVU)[5]='BioVU_odds_ratio'
names(BioVU)[6]='BioVU_P'
names(BioVU)[7]='BioVU_L95'
names(BioVU)[8]='BioVU_U95'
names(BioVU)[26]='BioVU_Power'

d=merge(d,BioVU[,c('assoc_ID','BioVU_cases','BioVU_controls','BioVU_odds_ratio','BioVU_P','BioVU_L95','BioVU_Power')], by='assoc_ID', all.x=T)
head(d)

table(d$phecode_string)


