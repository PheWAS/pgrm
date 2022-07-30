library(pgrm)
library(data.table)
library(glue)

results_annotated

head(results_annotated)

biovu_EUR=annotate_results(results_BioVU_EUR, build="hg19",ancestry="EUR")
get_RR(biovu_EUR) # Replicated 643 of 833 for RR=77.2%
get_AER(biovu_EUR) # Expected 1656.8, replicated 1354 for AE=0.817 (3268 associations for 106 uniq phecodes)
table(biovu_EUR$CI_overlap)
