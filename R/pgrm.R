#' @import checkmate
#' @import data.table
#' @import gaston
#' @import genpwr
#' @importFrom glue glue
NULL


#' Get instance of PGRM
#'
#' This function generates a PGRM copy with specified ancestry, build, and phecode version
#'
#' @param ancestry A string that indicates the ancestry of the PGRM.
#' Options EAS, EUR, AFR, SAS, AMR, ALL. Default ALL
#' @param build A string indicating the genome reference build. Options hg19, hg38. Default is hg19.
#' @param phecode_version A string indicating the phecode version. Currently only
#' V1.2 is supported, which is the default
#' @param unique If TRUE, then rows are uniqued by SNP/phecode (only relevant when ancestry == 'ALL')
#'
#' @return A data.table of the PGRM
#'
#' @details This function returns a copy of the PGRM. The function assigns a column 'SNP' as
#' either SNP_hg19 or SNP_hg38, depending on the build argument. The returned PGRM only includes
#' associations annotated with the specified ancestry, and the risk allele frequency (RAF) is based
#' on the specified ancestry.
#'
#' @seealso [PGRM_ALL], [annotate_results()]
#'
#' @eval example_get_PGRM()
#'
#' @export
get_PGRM = function(ancestry = 'all', build = 'hg19', phecode_version = 'V1.2', unique = TRUE) {
  ## Avoid warnings about global vars
  cat_LOG10_P = SNP = . = phecode = NULL

  ancestry = toupper(ancestry)
  build = tolower(build)
  checkBuild(build)
  checkAncestry(ancestry)
  checkPhecodeVersion(phecode_version)

  PGRM_new = copy(pgrm::PGRM_ALL)
  if (build == 'hg19') {
    PGRM_new$SNP_hg38 = NULL
    setnames(PGRM_new, 'SNP_hg19', 'SNP')
  } else if (build == 'hg38') {
    PGRM_new$SNP_hg19 = NULL
    setnames(PGRM_new, 'SNP_hg38', 'SNP')
    }

  if (ancestry != 'ALL') {
    a = ancestry
    PGRM_new = PGRM_new[ancestry == a]
  } else if (ancestry == 'ALL' && isTRUE(unique)) {
    ## make the 'ALL' PGRM unique by SNP/phecode
    uniq_PGRM = PGRM_new[, .(cat_LOG10_P = max(cat_LOG10_P)), by = c('SNP', 'phecode')]
    PGRM_new = merge(uniq_PGRM, PGRM_new, by = c('SNP', 'phecode', 'cat_LOG10_P'))}

  freq_col_name = glue('{ancestry}_RAF')
  setnames(PGRM_new, freq_col_name, 'RAF')
  PGRM_new = PGRM_new[, c(
   'assoc_ID', 'SNP', 'ancestry', 'rsID', 'risk_allele_dir', 'RAF',
   'phecode', 'phecode_string', 'category_string', 'cat_LOG10_P', 'cat_OR', 'cat_L95',
   'cat_U95', 'Study_accession', 'pub_count', 'pub_date')]

  return(PGRM_new)}


#' Annotate a result set with the PGRM
#'
#' This function annotates a result from a test cohort with information from the PGRM
#'
#' @param results A data frame (or data.table) with results of a test cohort; columns
#' for SNP, phecode, cases, controls, odds_ratio, P (see demo files for example (e.g results_MGI))
#' @param use_allele_dir If TRUE, direction of effect is used when assessing if an
#' association is replicated. To use this argument, odds ratios must be reported
#' for the alternative allele
#' @param ancestry A string that specifies ancestry of the PGRM that is then used
#' to annotate the results file.
#' Options EAS, EUR, AFR, SAS, AMR, ALL. Default ALL
#' @param build A string indicating the genome reference build used in the results
#' table. Options hg19, hg38. Default is hg19.
#' @param phecode_version A string indicating the phecode version used in the results
#' table. (Currently only V1.2 is supported)
#' @param calculate_power If TRUE then power calculations will be conducted using
#' case and control counts from the results file. Necessary for get_AER(). Default FALSE
#' @param annotate_CI_overlap If TRUE then a column called 'annotate_CI_overlap'
#' is added to the table, values:
#' **overlap**: 95% CIs of PGRM and test cohort overlap
#' **test_cohort_greater**: 95% CI of test cohort greater than PGRM
#' **PGRM_greater**: 95% CI of PGRM greater than test cohort
#' If annotate_CI_overlap is TRUE, then results must include 95% CIs
#' @param LOUD If TRUE then progress info is printed to the terminal. Default FALSE
#'
#' @return A data.table of the results file annotated with columsn from the PGRM
#'
#' @details This function takes a dataframe with summary statistics from a test cohort.
#' For an example of the way to format the results data frame, see one of the results
#' sets included in the package (e.g. results_MGI). (NOTE: If the direction of effect
#' is used to determine if an association is replicated, then the odds ratios of
#' the result set must be oriented to the alternative allele.)
#'
#' The function returns a data.table with the following annotations:
#' \itemize{
#'   \item Phecode informtion, including phecode_string and phecode_category
#'   \item Risk allele frequency from GnomAD (column RAF), ancestry specified by the ancestry argument
#'   \item The rsID
#'   \item The direction of effect (ref or alt) and risk allele of the original association
#'   \item Summary statistics from the GWAS catalog association, including the -log10(P),
#'   odds ratio, and 95% confidence intervals (cat_LOG10_P, cat_OR, cat_L95, cat_U95)
#'   \item The study accession ID from the GWAS catalog
#'   \item A column called powered, which is 1 or 0 indicating whether the test association
#'   is powered > 80%. If calculate_power==TRUE, then the power is determined by the case/control
#'   counts specified in the results data.table. Otherwise, it is derrived from the estimate
#'   pre-computed cases needed assuming a 1:5 case:control ratio. All power calculations use
#'   alpha=0.0
#'   \item A column called rep that indicates if the association is replicated (i.e.
#'   p<0.05 in the test cohort; if use_allele_dir==TRUE, then the direction of effect
#'   from the test cohort must also be consistant with what is reported in the catalog)
#'   \item If annotate_CI_overlap is true, then information about the relationship
#'   between the 95% CIs from the catalog and the test set is included in column
#'   CI_overlap, and new columns for odds_ratio, L95, and U95 are created(rOR, rL95,
#'   rU95) that are oriented to the risk allele. (This option requires that the confidence
#'   intervals are reported in the test cohort summary statistics)
#' }
#'
#' @eval annotate_UKBB()
#'
#' @seealso Example result sets: [results_BBJ], [results_UKBB], [results_BioVU_EUR], [results_BioVU_AFR], [results_MGI]
#'
#' @export
annotate_results = function(
    results, use_allele_dir = TRUE, ancestry = 'all', build = 'hg19', phecode_version = 'V1.2',
    calculate_power = TRUE, annotate_CI_overlap = TRUE, LOUD = TRUE) {
  ## Avoid warnings about global vars
  cases = cases_needed = risk_allele_dir = odds_ratio = cat_L95 = cat_U95 = rL95 =
    rU95 = powered = P = rOR = L95 = U95 = CI_overlap = SNP = Power = NULL
  PGRM = get_PGRM(ancestry = ancestry, build = build, phecode_version = phecode_version)
  checkResults(results)
  if (isTRUE(annotate_CI_overlap)) {
    checkForCIs(results)}

  results = data.table(results)
  results[, SNP := toupper(SNP)]
  results = merge(results, PGRM, by = c('SNP', 'phecode'))

  results[, powered := 0]
  results[!is.na(cases_needed) & cases >= cases_needed, powered := 1]

  results[, rep := 0]
  results[P < 0.05, rep := 1]

  if (use_allele_dir) {
    results[risk_allele_dir == 'ref' & odds_ratio > 1, rep := 0]
    results[risk_allele_dir == 'alt' & odds_ratio < 1, rep := 0]}
  if (calculate_power == TRUE) {
    results = annotate_power(results, LOUD = LOUD)
    results[, powered := 0]
    results[!is.na(Power) & Power >= 0.8, powered := 1]}
  if (annotate_CI_overlap == TRUE) {
    results[, rOR := odds_ratio]
    results[, rL95 := L95]
    results[, rU95 := U95]
    results[risk_allele_dir == 'ref', rOR := 1 / odds_ratio]
    results[risk_allele_dir == 'ref', rU95 := 1 / L95]
    results[risk_allele_dir == 'ref', rL95 := 1 / U95]

    results[, CI_overlap := '']
    results[rL95 >= cat_L95 & rL95 <= cat_U95, CI_overlap := 'overlap']
    results[rU95 >= cat_L95 & rL95 <= cat_U95, CI_overlap := 'overlap']
    results[rL95 > cat_U95, CI_overlap := 'test_cohort_greater']
    results[cat_L95 > rU95, CI_overlap := 'PGRM_greater']}

  setcolorder(results, results[, c('assoc_ID', names(results)[names(results) != 'assoc_ID'])])

  return(results)}


#' Calculate the replicaiton rate (RR) of a test cohort
#'
#' This function calculates the replicaiton rate in a test cohort that has been annotated with PGRM.
#' By default, it calculates the RR for associations that powered at >80%.
#'
#' @param annotated_results A data table of results that have been annotated with the PGRM
#' @param include A character string. If 'powered' then only powered associations are included.
#' If 'all' then all associations are included
#' @param LOUD If TRUE then progress info is printed to the terminal. Default TRUE
#'
#' @return An numeric value of the replication rate of the result set
#'
#'
#' @export
get_RR = function(annotated_results, include = 'powered', LOUD = TRUE) {
  ## Avoid warnings about global vars
  powered = NULL

  include = tolower(include)
  check_RR_include(include)
  checkAnnotatedResults(annotated_results, include)

  if (include == 'powered') {
    denominator = nrow(annotated_results[powered == 1])
    numerator = nrow(annotated_results[powered == 1 & rep == 1])
  } else {
    denominator = nrow(annotated_results)
    numerator = nrow(annotated_results[rep == 1])}
  RR = numerator / denominator
  if (LOUD == TRUE) {
    print(glue('Replicated {numerator} of {denominator} for RR={f_RR}',
               f_RR = sprintf('%1.1f%%', 100 * RR)))}
  return(RR)}


#' Calculate the acual:expected ratio (A:E) of a test cohort
#'
#' This function calculates the actual:expected a test cohort that has been annotated with PGRM.
#'
#' @param annotated_results A data table of results that have been annotated with the PGRM
#' @param LOUD If TRUE then progress info is printed to the terminal. Default TRUE
#'
#' @return An numeric value of the actual:expected ratio of the result set
#'
#' @export
get_AER = function(annotated_results, LOUD = TRUE) {
  ## Avoid warnings about global vars
  Power = NULL

  checkAnnotatedResults_forAE(annotated_results)

  r = annotated_results[!is.na(Power)]
  expected = sum(r$Power)
  actual = sum(r$rep)
  AE = actual / expected
  expected = round(expected, 1)
  total_assoc = nrow(r)
  uniq_phecode = length(unique(r$phecode))
  if (LOUD == TRUE) {
    ## Need to keep this line long because it prints to the terminal
    print(glue('Expected {expected}, replicated {actual} for AE={AE_round} ({total_assoc} associations for {uniq_phecode} uniq phecodes)',
               AE_round = round(AE, 3)))}
  return(AE)}


#' Calculate the percent of associations that are powered in a test set
#'
#' This function calculates % of powered associations for a test cohort that has been annotated with PGRM.
#'
#' @param annotated_results A data table of results that have been annotated with the PGRM
#' @param include_missing_pheno If TRUE, then denominator inclues SNP/phecode pairs that are not
#' present in the annotated_results data.table
#' @param LOUD If TRUE then progress info is printed to the terminal. Default TRUE
#'
#' @return An numeric value of the percent of associations that are powered
#'
#' @export
get_powered_rate = function(annotated_results, include_missing_pheno = TRUE, LOUD = TRUE) {
  checkAnnotatedResults(annotated_results, include = 'powered')

  total_rows = nrow(annotated_results)
  powered = nrow(annotated_results[powered == 1])
  powered_rate = powered / total_rows
  pr = round(powered_rate, 2)

  if (LOUD == TRUE) {
    uniq_phecode = length(unique(annotated_results[powered == 1]$phecode))
    ## Need to keep this line long because it prints to the terminal
    print(glue('Powered for {powered} of {total_rows} associations {pr} for {uniq_phecode} uniq phecodes'))}
  return(powered_rate)}
