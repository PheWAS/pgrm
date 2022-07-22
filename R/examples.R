example_get_PGRM = function() {
    ex = "
@examples
library(pgrm)

## Get a copy of the PGRM for build hg19, East Asian ancestry
get_PGRM(build='hg19',ancestry='EAS')
"
  return(strsplit(ex, split = '\n')[[1L]])}

annotate_UKBB = function() {
  ex = "
@examples
library(pgrm)

## annotate the BioVU African ancestry result set
anno = annotate_results(results_BioVU_AFR,ancestry='AFR', build='hg19',calculate_power=TRUE)

## Get the replication rate of associations powered at >80%
get_RR(anno)

## Get the replication rate of all associations
get_RR(anno,include='all')

## Get the actual:expected ratio
get_AER(anno)
"
  return(strsplit(ex, split = '\n')[[1L]])}
