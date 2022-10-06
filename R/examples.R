example_get_PGRM = function() {
    ex = "
@examples
library(pgrm)

## Get a copy of the PGRM for build hg19, East Asian ancestry
get_PGRM(build = 'hg19', ancestry = 'eas')
"
  return(strsplit(ex, split = '\n')[[1L]])}

annotate_UKBB = function() {
  ex = "
@examples
library(pgrm)

## get summary statistics for UK Biobanks
summary_stats_ukb = get_summary_stats(dataset = 'ukb')

## Filter out associations where UKB was included in the original study
summary_stats_ukb = summary_stats_ukb[cohort_match == 0]

## annotate summary statistics with pgrm
anno = annotate_results(head(summary_stats_ukb,n=10), ancestry = 'EUR',
  build = 'hg38', calculate_power = TRUE)

## Get the replication rate of associations powered at >80%
get_RR(anno)

## Get the replication rate of all associations
get_RR(anno,include='all')

## Get the actual:expected ratio
get_AER(anno)
"
  return(strsplit(ex, split = '\n')[[1L]])}

#' run_PGRM_assoc_ex = function (){
#'   ex = "
#' @examples
#'   library(pgrm)
#'   library(data.table)
#'
#'   ## Read in an ICD file. Format should look like the example file.
#'   ## Note: Don't forget to load icd column as a character vector!
#'   head(icdExampleTable)
#'
#'   ## make pheno table
#'   pheno = make_pheno(icdExampleTable)
#'   head(pheno)
#'
#'   ## Load in covariate file. Here's a generated example
#'   covar_table_test = data.table(
#'   person_id = 1:4, sex = c('F', 'M', 'F', 'M'),
#'   last_age = c(28772, 18028, 11636, 14589))
#'
#'   ## Load genotype file
#'   #library(gaston)
#'   #geno = read.bed.matrix('data/geno_test')
#'
#'   ## Make a PGRM instance, specifying genome build and ancestry
#'   #PGRM = get_PGRM(build = 'hg19', ancestry = 'all')
#'
#'   ## Run associations from the PGRM
#'  # run_PGRM_assoc(geno, pheno, covars, covariate_names = c('last_age'),
#' # PGRM, MCC = 2, use_exclude_range = TRUE, check_sex = TRUE)
#'   "
#'   return(strsplit(ex, split = '\n')[[1L]])}
