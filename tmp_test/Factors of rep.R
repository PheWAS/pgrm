library(PGRM)
library(splines)

d=benchmark_results

d$ancestry = as.factor(d$ancestry)
d$ancestry = relevel(d$ancestry, ref="EUR")

d$category_string = as.factor(d$category_string)
d$category_string = relevel(d$category_string, ref="neoplasms")

d$cohort = as.factor(d$cohort)
d$cohort = relevel(d$cohort, ref="MGI")






m=glm(rep~cases+cohort, data=d,family="binomial")
summary(m)
confint.default(m)
plot(allEffects(m))
tab_model(m)

d$logP = -1*log10(d$P)

m=glm(logP~Power+category_string+ancestry+first_pub, data=d)
summary(m)
plot_model(m)
tab_model(m)

m=glm(rep~cases+log(cat_L95)+category_string+ancestry+first_pub, data=d,family="binomial")
summary(m)
plot_model(m)
tab_model(m)



d$cat_L95_trans = log(d$cat_L95)*10
m2=glm(rep~cases+cat_L95_trans+AF+ancestry, data=d,family="binomial")
summary(m2)


m2=glm(rep~cat_L95_trans+cases+AF+ancestry+category_string+first_pub, data=d,family="binomial")
#m2=glm(rep~category_string+ancestry, data=d,family="binomial")
summary(m2)
exp(coef(m2)["cases"]*1000)
exp(coef(m2)["first_pub"]*365)
plot_model(m2)
confint.default(m2)
tab_model(m2)

m4=glm(P~cases+cat_LOG10_P+cat_L95+category_string+ancestry, data=d)
summary(m4)
plot_model(m4)
plot(allEffects(m))
library(sjPlot)
library(sjmisc)
library(sjlabelled)

tab_model(m4)

biovu_EUR$BioVU = 1
biovu_EUR$BioVU_rep = 0
biovu_EUR[rep==1]$BioVU_rep = 1

anno_MGI$MGI = 1
anno_MGI$MGI_rep = 0
anno_MGI[rep==1]$MGI_rep = 1

anno_UKBB$UKBB=1
anno_UKBB$UKBB_rep=0
anno_UKBB[rep==1]$UKBB_rep = 1


#d=PGRM_ALL[,c('assoc_ID','phecode','phecode_string','category_string','ancestry','cat_LOG10_P','cat_OR', 'cat_L95','cat_U95')]
#head(biovu_EUR)


#d=merge(biovu_EUR[powered==1,c("assoc_ID","BioVU","BioVU_rep")],anno_MGI[powered==1,c("assoc_ID","MGI","MGI_rep")],by="assoc_ID")
#d=merge(d,anno_UKBB[powered==1,c("assoc_ID","UKBB","UKBB_rep")],by="assoc_ID")
#nrow(d) ## 379
#head(d)


## Venn diagram

x <- list (
  BioVU =  d[d$BioVU_rep==1,]$assoc_ID,
  MGI =  d[d$MGI_rep==1,]$assoc_ID,
  UKBB =  d[d$UKBB_rep==1,]$assoc_ID
  #,null =  d[d$null==1,]$assoc_ID
)
library(ggvenn)
ggvenn(
  x,
  fill_color = c("#EFC000FF", "#0073C2FF", "#CD534CFF", "#CD534CFF"),
  stroke_size = 0.8, set_name_size = 8, show_percentage=TRUE,digits=0,text_size=0
)

table(PGRM_ALL[assoc_ID %in% d$assoc_ID]$category_string)
table(PGRM_ALL[assoc_ID %in% d[BioVU_rep==0 & MGI_rep==0 & UKBB_rep == 0]$assoc_ID]$category)
#table(PGRM_ALL[assoc_ID %in% d[BioVU_rep==0 & MGI_rep==0 & UKBB_rep == 0]$assoc_ID]$phecode_string)
#table(PGRM_ALL$category_string)

#psych==5

table(PGRM_ALL[assoc_ID %in% d$assoc_ID]$category)


## by category

biovu_EUR[, .(power_rep=sum(d$rep[d$powered==1]),Powered=sum(Power),Expected=sum(Power),actual=sum(rep)), by = "category_string"]

biovu_cat=biovu_EUR[, .(power_rep=uniqueN(assoc_ID[powered==1 & rep==1]),powered=sum(powered),RR_power=uniqueN(assoc_ID[powered==1 & rep==1])/sum(powered),
                        total_rep=sum(rep),Expected=sum(Power),AER=sum(rep)/sum(Power)), by = "category_string"]

MGI_cat=anno_MGI[, .(power_rep=uniqueN(assoc_ID[powered==1 & rep==1]),powered=sum(powered),RR_power=uniqueN(assoc_ID[powered==1 & rep==1])/sum(powered),
                     total_rep=sum(rep),Expected=sum(Power),AER=sum(rep)/sum(Power)), by = "category_string"]

UKBB_cat=anno_UKBB[!is.na(Power), .(power_rep=uniqueN(assoc_ID[powered==1 & rep==1]),powered=sum(powered),RR_power=uniqueN(assoc_ID[powered==1 & rep==1])/sum(powered),
                                    total_rep=sum(rep),Expected=sum(Power),AER=sum(rep)/sum(Power)), by = "category_string"]
UKBB_cat


biovu_afr_cat=biovu_AFR[, .(power_rep=uniqueN(assoc_ID[powered==1 & rep==1]),powered=sum(powered),RR_power=uniqueN(assoc_ID[powered==1 & rep==1])/sum(powered),
                            total_rep=sum(rep),Expected=sum(Power),AER=sum(rep)/sum(Power)), by = "category_string"]

BBJ_cat=anno_BBJ[, .(power_rep=uniqueN(assoc_ID[powered==1 & rep==1]),powered=sum(powered),RR_power=uniqueN(assoc_ID[powered==1 & rep==1])/sum(powered),
                     total_rep=sum(rep),Expected=sum(Power),AER=sum(rep)/sum(Power)), by = "category_string"]

biovu_cat=as.data.frame.matrix(table(biovu_EUR[powered==1]$category_string, biovu_EUR[powered==1]$rep))
#biovu_cat$BioVU_rep = biovu_cat$`1`/(biovu_cat$`1`+biovu_cat$`0`)
biovu_cat$cat = rownames(biovu_cat)
biovu_cat

mgi_cat=as.data.frame.matrix(table(anno_MGI[powered==1]$category_string, anno_MGI[powered==1]$rep))
mgi_cat$MGI_rep = mgi_cat$`1`/(mgi_cat$`1`+mgi_cat$`0`)
mgi_cat$cat = rownames(mgi_cat)
mgi_cat

ukbb_cat=as.data.frame.matrix(table(anno_UKBB[powered==1]$category_string, anno_UKBB[powered==1]$rep))
ukbb_cat$UKBB_rep = ukbb_cat$`1`/(ukbb_cat$`1`+ukbb_cat$`0`)
ukbb_cat$cat= rownames(ukbb_cat)

d=merge(biovu_cat, mgi_cat,by="cat")
d=merge(d,ukbb_cat,by="cat")
d

?pheatmap
library(pheatmap)

foo=as.matrix(d[,c("BioVU_rep","MGI_rep","UKBB_rep")])
rownames(foo)=d$cat
pheatmap(foo,display_numbers=T,cluster_row=F,cluster_cols = F)

test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")

