library(pgrm)

d=benchmark_results
table(d$cohort)
## by category

table(d[cohort=="BioVU_EUR"]$powered)
table(d[cohort=="BioVU_EUR" & powered == 1]$rep)
sum(d[cohort=="BioVU_EUR"]$rep)

biovu_cat=d[cohort=="BioVU_EUR", .(power_rep=uniqueN(assoc_ID[powered==1 & rep==1]),powered=sum(powered),RR_power=uniqueN(assoc_ID[powered==1 & rep==1])/sum(powered),
                        total_rep=sum(rep),Expected=sum(Power),AER=sum(rep)/sum(Power)), by = "category_string"]
biovu_cat[order(category_string)]

MGI_cat=d[cohort=="MGI", .(power_rep=uniqueN(assoc_ID[powered==1 & rep==1]),powered=sum(powered),RR_power=uniqueN(assoc_ID[powered==1 & rep==1])/sum(powered),
                     total_rep=sum(rep),Expected=sum(Power),AER=sum(rep)/sum(Power)), by = "category_string"]

MGI_cat[order(category_string)]

UKBB_cat=d[cohort=="UKBB", .(power_rep=uniqueN(assoc_ID[powered==1 & rep==1]),powered=sum(powered),RR_power=uniqueN(assoc_ID[powered==1 & rep==1])/sum(powered),
                                    total_rep=sum(rep),Expected=sum(Power),AER=sum(rep)/sum(Power)), by = "category_string"]
UKBB_cat[order(category_string)]


biovu_afr_cat=d[cohort=="BioVU_AFR", .(power_rep=uniqueN(assoc_ID[powered==1 & rep==1]),powered=sum(powered),RR_power=uniqueN(assoc_ID[powered==1 & rep==1])/sum(powered),
                            total_rep=sum(rep),Expected=sum(Power),AER=sum(rep)/sum(Power)), by = "category_string"]

biovu_afr_cat[order(category_string)]

BBJ_cat=d[cohort=="BBJ", .(power_rep=uniqueN(assoc_ID[powered==1 & rep==1]),powered=sum(powered),RR_power=uniqueN(assoc_ID[powered==1 & rep==1])/sum(powered),
                     total_rep=sum(rep),Expected=sum(Power),AER=sum(rep)/sum(Power)), by = "category_string"]
BBJ_cat[order(category_string)]

