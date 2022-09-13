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
#' @eval run_PGRM_assoc_ex()
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
    print(glue('Expected {expected}, replicated {actual} for AE={AE_round} ({total_assoc}
               associations for {uniq_phecode} uniq phecodes)', AE_round = round(AE, 3)))}
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
    print(glue('Powered for {powered} of {total_rows} associations {pr} for {uniq_phecode}
               uniq phecodes'))}
  return(powered_rate)}


#' Make a phecode phenotype file from table of ICD codes
#'
#' This function calculates % of powered associations for a test cohort that has been annotated with PGRM.
#'
#' @param icds A data table ICD codes. Must include columns `person_id`, `icd`, `flag`
#' @param phecode_version Phecode map version. Currently only supports 'V1.2', which is default
#'
#' @return An numeric value of the percent of associations that are powered
#'
#' @export
make_pheno = function(icds, phecode_version = 'V1.2') {
  checkIcdTable(icds)
  checkPhecodeVersion(phecode_version)
  if(phecode_version == 'V1.2') {
    icdPhecodeMap = icdPhecodeMap_V1_2}
  pheno = merge(icds, icdPhecodeMap, by = c('icd', 'flag'), allow.cartesian = TRUE)
  pheno = pheno[, c('person_id', 'phecode', 'entry_date')]
  pheno = unique(pheno)
  pheno = pheno[, .N, by = c('person_id', 'phecode')]
  return(pheno)}

#' Add case/control status for a specified phenotype to a covariate file.
#'
#' @param pheno A data.table with phecode phenotypes. columns `person_id`, `phecode`, and `N`
#' required. The phecode columns should be a character vector. Column `N` specifies the number of
#' unique dates a phecode occurs for that individual. Only phecodes that are present for an individual
#' are included (i.e. N is always greater than one).
#' @param covars A data.table having one row per person in the cohort. Must havea column names
#' `person_id`. Must include other columns specified in covariate_names.
#' @param phecode A character vector of the phecode
#' @param MCC A numeric value specified the minimum code count required for cases.
#' @param use_exclude_ranges If TRUE then individuals who are not cases but have a phecode within
#' the exclude ranges are excluded.
#' @param check_sex If TRUE then the `pheno` column is set to zero for individuals with a sex that
#' does not match the specified sex. This argument is only relevant for sex specific codes
#' (e.g. Prostate cancer). To use this function, sex must be included as a column in the covars data.table.
#'
#' @return The covar data.table with a new column called `pheno`. This column is 0 if the individual,
#' is a control for specified phecode, 1 if they are a case, and NA if they are excluded.
#'
#' @export
get_pheno = function(pheno, covars, phecode, MCC = 2, use_exclude_ranges = TRUE, check_sex = FALSE) {
  checkPhecodeTable(pheno)
  checkDemosTable(covars)
  checkPhecode(phecode)
  checkSex(check_sex, covars)
  checkMCC(MCC)

  phecode_num = N = sex = NULL

  cur_phecode = phecode

  p = pheno[pheno$phecode == cur_phecode, c('person_id', 'N')]
  p = data.table(p, key = 'person_id')
  p$pheno = -9
  p[p$N >= MCC]$pheno = 1

  if(use_exclude_ranges == TRUE){
    to_exclude = c()
    pheno$phecode_num = as.numeric(pheno$phecode)
    ranges = pgrm::exclude_ranges[phecode == cur_phecode]
    for(i in 1:nrow(ranges)){
      to_exclude = union(to_exclude, unique(pheno[phecode_num >= ranges[i,]$range_start &
                                                    phecode_num < ranges[i,]$range_end]$person_id))
    }
    to_exclude = subset(to_exclude, !to_exclude %in% p$person_id)
    to_exclude = data.table(person_person_id = to_exclude)
    names(to_exclude)[1] = 'person_id'
    to_exclude$N = 0
    to_exclude$pheno = -9
    p = rbind(p, to_exclude)
  }
  p = merge(covars, p, by = 'person_id', all.x=TRUE)
  p[is.na(pheno)]$pheno = 0
  p[p$pheno == -9]$pheno = NA
  p[is.na(N)]$N = 0

  if(check_sex == TRUE){
    phecode_sex = sex_check_phecode(cur_phecode)
    if(phecode_sex %in% c('F', 'M')){
      p[sex != phecode_sex]$pheno = NA
    }
  }
  return(p)}

#' Run associations in the pgrm on data from a test cohort. Function requires row level phenotype,
#' genotype, and covariate information.
#'
#' (NOTE: If you used principle components as covariates, make sure they are not coded as numeric.
#' Otherwise the logistic regression will treat them as categorical variables, making the function
#' run very slowly before failing.)
#'
#' @param geno A matrix or 'BEDMatrix' object containing genetic data
#' @param pheno A data.table with phecode phenotypes. columns `person_id`, `phecode`, and `N`
#' required. The phecode columns should be a character vector. Column `N` specifies the number of
#' unique dates a phecode occurs for that individual. Only phecodes that are present for an individual
#' are included (i.e. N is always greater than one).
#' @param covars A data.table having one row per person in the cohort. Must havea column names
#' `person_id`. Must include other columns specified in covariate_names.
#' @param covariate_names A list of covariates to be used in the logistic regression. All covariates
#'   listed must be present in the covars table
#' @param PGRM A data.table of the PGRM, generated with get_PGRM()
#' @param MCC An integer specifying the minimum code count required for a case (default MCC=2)
#' @param minimum_case_count An integer specifiying the minimum number of cases required in the test
#'   cohort to be included in the analysis (default minimum_case_count=100)
#' @param use_exclude_ranges If TRUE then exclude ranges are applied to controls
#' @param check_sex If TRUE then individuals with sex == F in the covars table will be excluded
#'   from male-specific phenotypes, and vice-versa. Sex specific phecodes are specified in
#'   phecode_info. For this function to work, the covars table must have column `sex` and the values
#'   must be 'M' for Male and 'F' for Female
#' @param LOUD If TRUE then progress info is printed to the terminal. Default TRUE
#'
#' @return A data.table with annotated results from association tests
#'
#' @eval run_PGRM_assoc_ex()
#'
#' @export
run_PGRM_assoc = function(geno, pheno, covars,covariate_names, PGRM, MCC=2, minimum_case_count = 100,
                          use_exclude_ranges = TRUE, check_sex = FALSE, LOUD = TRUE) {
  person_id = id = N = . = phecode = cases = NULL
  ## Add check PGRM
  checkGenotypes(geno)
  checkPhecodeTable(pheno)
  checkDemosTable(covars)
  checkCovarList(covariate_names,covars)
  checkMCC(MCC)
  checkMCC(minimum_case_count)
  checkBool(use_exclude_ranges)
  checkSex(check_sex, covars)
  checkBool(LOUD)
  # ## create formulas for glm using covariates; formula_string_no_sex is for sex-specific phenotypes
  covar_list = paste(covariate_names, collapse = '+')
  formula_string = paste('pheno~genotype+', covar_list, sep = '')
  formula_string_no_sex = gsub('\\+sex', '', formula_string)
  formula_string_no_sex = gsub('sex\\+', '', formula_string_no_sex)

  ## Change person_id to character to ensure merging works
  covars$person_id = as.character(covars$person_id)
  pheno$person_id = as.character(pheno$person_id)

  # filter covars, pheno, and geno for intersection of IDs
  IDs = union(geno@ped$id, covars$person_id)
  covars = covars[person_id %in% IDs]
  geno = select.inds(geno, id %in% IDs)
  pheno = pheno[person_id %in% IDs]

  ## get available phecodes list
  phecode_counts = pheno[N >= MCC, .(cases = .N), by = phecode]
  available_phecodes = phecode_counts[cases >= minimum_case_count]$phecode

  # filter PGRM for available SNPs and phecodes
  SNPs = geno@snps$id
  PGRM = PGRM[PGRM$SNP %in% SNPs & PGRM$phecode %in% available_phecodes,]

  assoc_to_run = unique(PGRM[, c('SNP', 'phecode')])
  assoc_num = nrow(assoc_to_run)
  if(assoc_num == 0){
    print('There are no eligable associations to run')
    return(0)}
  if(LOUD == TRUE){
    print(glue('Attempting {assoc_num} association tests'))
    print(glue('Formula: {formula_string}'))}

  results = data.frame()

  for(i in 1:nrow(assoc_to_run)){
    cur_SNP = assoc_to_run[i,]$SNP
    cur_phecode = assoc_to_run[i,]$phecode
    cur_SNP_index = which(SNPs %in% c(cur_SNP))

    if(LOUD == TRUE){
      print(glue('{i} SNP: {cur_SNP} Phecode: {cur_phecode}'))
    }

    g = data.table(as.matrix(geno[, cur_SNP_index]), keep.rownames = TRUE)
    names(g)[1] = 'person_id'
    names(g)[2] = 'genotype'
    setkey(g, 'person_id')

    g$genotype = abs(g$genotype - 2) ## gaston codes 2==HOM for ref. Flip this around

    d = merge(g, covars, by='person_id')
    d = get_pheno(pheno = pheno, covars = d, phecode = cur_phecode, MCC = MCC,
                use_exclude_ranges = use_exclude_ranges, check_sex = check_sex)

    formula = as.formula(formula_string)
    if(check_sex == TRUE) {
      phecode_sex = sex_check_phecode(cur_phecode)
      if(phecode_sex != 'B') {
        formula = as.formula(formula_string_no_sex)}}

    n_case = nrow(d[pheno == 1,])
    n_control = nrow(d[pheno == 0,])
    if(n_case < minimum_case_count) {
      next}

    m = glm(formula, data = d, family = 'binomial')
    conf = confint.default(m)
    P = summary(m)$coeff[2, 4]
    odds_ratio = exp(summary(m)$coeff[2, 1])

    L95 = exp(conf[2, 1])
    U95 = exp(conf[2, 2])
    result = data.frame(SNP = cur_SNP, phecode = cur_phecode, cases = n_case, controls = n_control,
                        P = P, odds_ratio = odds_ratio, L95 = L95, U95 = U95)
    results = rbind(results, result)}
  return(results)}
