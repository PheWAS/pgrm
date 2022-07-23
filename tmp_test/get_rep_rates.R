library(pgrm)
library(data.table)


biovu_EUR=annotate_results(results_BioVU_EUR[cohort_match==0], build="hg19",ancestry="EUR")
get_RR(biovu_EUR)
get_AER(biovu_EUR)

biovu_AFR=annotate_results(results_BioVU_AFR[cohort_match==0], build="hg19",ancestry="AFR")
get_RR(biovu_AFR)
get_AER(biovu_AFR)

anno_BBJ=annotate_results(results_BBJ[cohort_match==0], build="hg19",ancestry="EAS")
get_RR(anno_BBJ)
get_AER(anno_BBJ)


anno_MGI=annotate_results(results_MGI[cohort_match==0], build="hg38",ancestry="EUR")
get_RR(anno_MGI)
get_AER(anno_MGI)

anno_UKBB=annotate_results(results_UKBB[cohort_match==0], build="hg38",ancestry="EUR")
get_RR(anno_UKBB)
get_AER(anno_UKBB)

biovu_EUR$BioVU = 1
biovu_EUR$BioVU_rep = 0
biovu_EUR[rep==1]$BioVU_rep = 1

anno_MGI$MGI = 1
anno_MGI$MGI_rep = 0
anno_MGI[rep==1]$MGI_rep = 1

anno_UKBB$UKBB=1
anno_UKBB$UKBB_rep=0
anno_UKBB[rep==1]$UKBB_rep = 1

d=merge(biovu_EUR[powered==1,c("assoc_ID","BioVU","BioVU_rep")],anno_MGI[powered==1,c("assoc_ID","MGI","MGI_rep")],by="assoc_ID")
d=merge(d,anno_UKBB[powered==1,c("assoc_ID","UKBB","UKBB_rep")],by="assoc_ID")
nrow(d)


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
  stroke_size = 0.8, set_name_size = 8, show_percentage=TRUE
)

table(PGRM_ALL[assoc_ID %in% d[BioVU_rep==0 & MGI_rep==0 & UKBB_rep == 0]$assoc_ID]$phecode_string)

table(PGRM_ALL[assoc_ID %in% d[BioVU_rep==0 & MGI_rep==0 & UKBB_rep == 0]$assoc_ID]$category)
#psych==5

table(PGRM_ALL[assoc_ID %in% d$assoc_ID]$category)


## by category

biovu_cat=as.data.frame.matrix(table(biovu_EUR[powered==1]$category_string, biovu_EUR[powered==1]$rep))
biovu_cat$BioVU_rep = biovu_cat$`1`/(biovu_cat$`1`+biovu_cat$`0`)
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

