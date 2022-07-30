library(pgrm)
library(data.table)
library(glue)
library(DescTools)

prep_table = function(anno, prefix=""){
  anno = anno[,c(10,3:9,24:26,30)]
  for(i in 2:ncol(anno)){
    names(anno)[i] = prefix %c% "_" %c% names(anno)[i]
  }
  return(anno)
}

biovu_EUR=annotate_results(results_BioVU_EUR, build="hg19",ancestry="EUR")
get_RR(biovu_EUR) # Replicated 643 of 833 for RR=77.2%
get_AER(biovu_EUR) # Expected 1656.8, replicated 1354 for AE=0.817 (3268 associations for 106 uniq phecodes)
table(biovu_EUR$CI_overlap)





biovu_AFR=annotate_results(results_BioVU_AFR, build="hg19",ancestry="AFR")
get_RR(biovu_AFR) ## Replicated 11 of 14 for RR=78.6%
get_AER(biovu_AFR) ## Expected 20, replicated 19 for AE=0.948 (31 associations for 14 uniq phecodes)

anno_BBJ=annotate_results(results_BBJ, build="hg19",ancestry="EAS")
get_RR(anno_BBJ) ## Replicated 167 of 218 for RR=76.6%
get_AER(anno_BBJ) ## Expected 259.6, replicated 205 for AE=0.79 (384 associations for 26 uniq phecodes)


anno_MGI=annotate_results(results_MGI, build="hg38",ancestry="EUR")
get_RR(anno_MGI) ## Replicated 694 of 905 for RR=76.7%
get_AER(anno_MGI) ## Expected 1907.5, replicated 1521 for AE=0.797 (4117 associations for 109 uniq phecodes)

anno_UKBB=annotate_results(results_UKBB, build="hg38",ancestry="EUR")
get_RR(anno_UKBB) # Replicated 706 of 818 for RR=86.3%
get_AER(anno_UKBB) # Expected 1327.5, replicated 1274 for AE=0.96 (2236 associations for 81 uniq phecodes)


PGRM=PGRM_ALL[,c('assoc_ID','Study_accession','SNP_hg19','SNP_hg38','risk_allele_dir','rsID','phecode','phecode_string','category_string','ancestry')]


d=merge(PGRM,prep_table(biovu_EUR,"BioVU_EUR"), by="assoc_ID",all.x=TRUE)
d=merge(d,prep_table(biovu_AFR,"BioVU_AFR"), by="assoc_ID",all.x=TRUE)
d=merge(d,prep_table(anno_MGI,"MGI"), by="assoc_ID",all.x=TRUE)
d=merge(d,prep_table(anno_UKBB,"UKBB"), by="assoc_ID",all.x=TRUE)
d=merge(d,prep_table(anno_BBJ,"BBJ"), by="assoc_ID",all.x=TRUE)
head(d)



write.table(d,file="all_results.txt",sep="\t",col.names = T,row.names = F)
