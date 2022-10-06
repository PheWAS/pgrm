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
#' Options afr, eas, eur, sas, amr or all. Default `all`
#' @param build A string indicating the genome reference build. Options hg19, hg38. Default is hg19.
#' @param phecode_version A string indicating the phecode version. Currently only
#' V1.2 is supported, which is the default
#' @param unique If TRUE, then rows are uniqued by SNP/phecode (only relevant when ancestry == 'ALL')
#'
#' @returns A data.table of the PGRM. For column details, see [pgrm_all].
#'
#' @details This function returns a copy of the PGRM. The function assigns a column 'SNP' as
#' either SNP_hg19 or SNP_hg38, depending on the build argument. The returned PGRM only includes
#' associations annotated with the specified ancestry, and the risk allele frequency (RAF) is based
#' on the specified ancestry.
#'
#' @seealso [pgrm_all], [annotate_results()]
#'
#' @eval example_get_PGRM()
#'
#' @export
get_PGRM = function(ancestry = 'all', build = 'hg19', phecode_version = 'V1.2', unique = TRUE) {
  ## Avoid warnings about global vars
  cat_log10_p = snp = . = phecode = NULL

  ancestry = tolower(ancestry)
  build = tolower(build)
  checkBuild(build)
  checkAncestry(ancestry)
  checkPhecodeVersion(phecode_version)

  PGRM_new = copy(pgrm::pgrm_all)
  if (build == 'hg19') {
    PGRM_new$snp_hg38 = NULL
    setnames(PGRM_new, 'snp_hg19', 'snp')
  } else if (build == 'hg38') {
    PGRM_new$snp_hg19 = NULL
    setnames(PGRM_new, 'snp_hg38', 'snp')
    }

  if (ancestry != 'all') {
    a = ancestry
    PGRM_new = PGRM_new[ancestry == a]
  } else if (ancestry == 'all' && isTRUE(unique)) {
    ## make the 'ALL' PGRM unique by SNP/phecode
    uniq_PGRM = PGRM_new[, .(cat_log10_p = max(cat_log10_p)), by = c('snp', 'phecode')]
    PGRM_new = merge(uniq_PGRM, PGRM_new, by = c('snp', 'phecode', 'cat_log10_p'))}

  freq_col_name = glue('{ancestry}_raf')
  setnames(PGRM_new, freq_col_name, 'raf')
  PGRM_new = PGRM_new[, c(
   'assoc_id', 'snp', 'ancestry', 'rsid', 'risk_allele_dir', 'raf',
   'phecode', 'phecode_string', 'category_string', 'cat_log10_p', 'cat_or', 'cat_l95',
   'cat_u95', 'study_accession', 'pub_count', 'pub_date')]

  return(PGRM_new)}

#' Retrieve summary statistics
#'
#' This function returns a data.table of summary statics from biobank cohorts
#'
#' @param dataset A string that indicates the biobank dataset
#' Options ukb, mgi, bbj, biovu_eur, biovu_afr. Default all
#' @param exclude_overlap If TRUE then associations that were originally sourced from the dataset
#' will be excluded
#'
#' @returns A data.table of the summary statistics
#'
#' @details This function retrieves summary statistics from biobank cohorts.
#'
#' @export
get_summary_stats = function(dataset = 'ukb', exclude_overlap = TRUE){
  ds = cohort_match = NULL
  ds = tolower(dataset)
  checkDataset(ds)
  checkBool(exclude_overlap)
  sum_stat = copy(pgrm::summary_stats)
  sum_stat = sum_stat[dataset == ds]
  if(exclude_overlap == TRUE){
    sum_stat = sum_stat[cohort_match == 0]
  }
  return(sum_stat)
}

#' Annotate a result set with the PGRM
#'
#' This function annotates a result from a test cohort with information from the PGRM
#'
#' @param results A data frame (or data.table) with results of a test cohort; columns
#' for SNP, phecode, cases, controls, odds_ratio, P.
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
#' @param annotate_ci_overlap If TRUE then a column called 'annotate_ci_overlap'
#' is added to the table, values:
#' **overlap**: 95% CIs of PGRM and test cohort overlap
#' **test_cohort_greater**: 95% CI of test cohort greater than PGRM
#' **PGRM_greater**: 95% CI of PGRM greater than test cohort
#' If annotate_ci_overlap is TRUE, then results must include 95% CIs
#' @param LOUD If TRUE then progress info is printed to the terminal. Default FALSE
#'
#' @return A data.table of the results file annotated with columsn from the PGRM
#'
#' @details This function takes a dataframe with summary statistics from a test cohort.
#' For an example of the way to format the results data frame, look at the fomatting from
#' [summary_stats]. (NOTE: If the direction of effect
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
#'   odds ratio, and 95% confidence intervals (cat_log10_p, cat_or, cat_l95, cat_u95)
#'   \item The study accession ID from the GWAS catalog
#'   \item A column called powered, which is 1 or 0 indicating whether the test association
#'   is powered > 80%. If calculate_power==TRUE, then the power is determined by the case/control
#'   counts specified in the results data.table. Otherwise, it is derrived from the estimate
#'   pre-computed cases needed assuming a 1:5 case:control ratio. All power calculations use
#'   alpha=0.0
#'   \item A column called rep that indicates if the association is replicated (i.e.
#'   p<0.05 in the test cohort; if use_allele_dir==TRUE, then the direction of effect
#'   from the test cohort must also be consistant with what is reported in the catalog)
#'   \item If annotate_ci_overlap is true, then information about the relationship
#'   between the 95% CIs from the catalog and the test set is included in column
#'   ci_overlap, and new columns for odds_ratio, l95, and u95 are created(risk_or, risk_l95,
#'   risk_u95) that are oriented to the risk allele. (This option requires that the confidence
#'   intervals are reported in the test cohort summary statistics)
#' }
#'
#' @eval annotate_UKBB()
#'
#' @seealso Example result sets: [summary_stats]
#'
#' @export
annotate_results = function(
    results, use_allele_dir = TRUE, ancestry = 'all', build = 'hg19', phecode_version = 'V1.2',
    calculate_power = TRUE, annotate_ci_overlap = TRUE, LOUD = TRUE) {
  ## Avoid warnings about global vars
  cases = cases_needed = risk_allele_dir = odds_ratio = cat_l95 = cat_u95 = risk_l95 =
    risk_u95 = powered = p = risk_or = l95 = u95 = ci_overlap = snp = power = NULL
  PGRM = get_PGRM(ancestry = ancestry, build = build, phecode_version = phecode_version)
  checkResults(results)
  if (isTRUE(annotate_ci_overlap)) {
    checkForCIs(results)}

  results = data.table(results)
  results[, snp := toupper(snp)]
  results = merge(results, PGRM, by = c('snp', 'phecode'))

  results[, powered := 0]

  results[, rep := 0]
  results[, powered := 0]
  results[, power := 0]
  results[, ci_overlap := '']
  results[, risk_or := 0]
  results[, risk_u95 := 0]
  results[, risk_l95 := 0]

  results[p < 0.05, rep := 1]


  if (use_allele_dir) {
    results[risk_allele_dir == 'ref' & odds_ratio > 1, rep := 0]
    results[risk_allele_dir == 'alt' & odds_ratio < 1, rep := 0]}
  if (calculate_power == TRUE) {
    results = annotate_power(results, LOUD = LOUD)
    results[, powered := 0]
    results[!is.na(power) & power >= 0.8, powered := 1]}
  if (annotate_ci_overlap == TRUE) {
    results[, risk_or := odds_ratio]
    results[, risk_l95 := l95]
    results[, risk_u95 := u95]
    results[risk_allele_dir == 'ref', risk_or := 1 / odds_ratio]
    results[risk_allele_dir == 'ref', risk_u95 := 1 / l95]
    results[risk_allele_dir == 'ref', risk_l95 := 1 / u95]

    results[, ci_overlap := '']
    results[risk_l95 >= cat_l95 & risk_l95 <= cat_u95, ci_overlap := 'overlap']
    results[risk_u95 >= cat_l95 & risk_l95 <= cat_u95, ci_overlap := 'overlap']
    results[risk_l95 > cat_u95, ci_overlap := 'test_cohort_greater']
    results[cat_l95 > risk_u95, ci_overlap := 'PGRM_greater']}

  #setcolorder(results, results[, c('assoc_id', names(results)[names(results) != 'assoc_id'])])
  setcolorder(results, c('assoc_id', 'snp', 'phecode', 'phecode_string', 'category_string', 'cases',
                         'controls', 'odds_ratio', 'p', 'l95', 'u95', 'power', 'powered', 'rep',
                         'ci_overlap', 'risk_or', 'risk_l95', 'risk_u95',
                         'study_accession', 'pub_count', 'pub_date', 'ancestry', 'rsid',
                         'risk_allele_dir', 'raf',
                         'cat_log10_p', 'cat_or', 'cat_l95', 'cat_u95' ))

  ## Get rid of power or ci_overlap columns if those were not computed
  if(calculate_power == FALSE){
    results$power = results$powered = NULL}

  if(annotate_ci_overlap == FALSE){
    results$ci_overlap = results$risk_l95 = results$risk_u95 = results$risk_or = NULL}

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


