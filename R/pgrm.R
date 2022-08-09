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
#' @param ancestry A string that indicates the ancestry of the PGRM. Options EAS, EUR, AFR, SAS, AMR, ALL. Default ALL
#' @param build A string indicating the genome reference build. Options hg19, hg38. Default is hg19.
#' @param phecode_version A string indicating the phecode version. Currently only
#' V1.2 is supported, which is the default
#' @param unique If TRUE, then rows are uniqued by SNP/phecode (only relevant when ancestry == "ALL")
#'
#' @return A data.table of the PGRM
#'
#' @details This function returns a copy of the PGRM with specified columns for SNP, AF, and cases_needed
#' according to the arguments specified. The function assigns a column "SNP" as either SNP_hg19 or SNP_hg38, depending
#' on the build argument. It also assigns the allele frequency (AF) and cases_needed columns according to the ancestry.
#' If the ancestry is != "ALL", then the funciton also subsets PGRM_ALL to include only associations that were based in
#' cohorts with the specified ancestry. This funciton is called inside annotate_results(),
#' which can be used to annoate the results of a specific test cohort.
#'
#' (NOTE: If you used principle components as covariates, make sure they are not coded as numeric.
#' Otherwise the logistic regression will treat them as categorical variables, making the function
#' run very slowly before failing.)
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
    ## make the "ALL" PGRM unique by SNP/phecode
    uniq_PGRM = PGRM_new[, .(cat_LOG10_P = max(cat_LOG10_P)), by = c('SNP', 'phecode')]
    PGRM_new = merge(uniq_PGRM, PGRM_new, by = c('SNP', 'phecode', 'cat_LOG10_P'))}

  freq_col_name = glue('{ancestry}_RAF')
  cases_needed_col_name = glue('cases_needed_{ancestry}')
  setnames(PGRM_new, freq_col_name, 'RAF')
  setnames(PGRM_new, cases_needed_col_name, 'cases_needed')
  PGRM_new = PGRM_new[, c(
   'assoc_ID', 'SNP', 'ancestry', 'rsID', 'risk_allele_dir', 'RAF',
   'phecode', 'phecode_string', 'category_string', 'cat_LOG10_P', 'cat_OR', 'cat_L95',
   'cat_U95', 'cases_needed', 'Study_accession','pub_count','first_pub_date')]

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
#' @param annotate_CI_overlap If TRUE then a column called "annotate_CI_overlap"
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

  setcolorder(results,results[ , c("assoc_ID",   names(results)[names(results) != "assoc_ID"])])

  return(results)}


#' Calculate the replicaiton rate (RR) of a test cohort
#'
#' This function calculates the replicaiton rate in a test cohort that has been annotated with PGRM.
#' By default, it calculates the RR for associations that powered at >80%.
#'
#' @param annotated_results A data table of results that have been annotated with the PGRM
#' @param include A character string. If "powered" then only powered associations are included.
#' If "all" then all associations are included
#' @param LOUD If TRUE then progress info is printed to the terminal. Default TRUE
#'
#' @return An numeric value of the replication rate of the result set
#'
#' @eval annotate_UKBB()
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
    print(glue('Expected {expected}, replicated {actual} for AE={AE_round} ({total_assoc} associations for {uniq_phecode} uniq phecodes)', AE_round = round(AE, 3)))}
  return(AE)}


#' Calculate the percent of associations that are powered in a test set
#'
#' This function calculates % of powered associations for a test cohort that has been annotated with PGRM.
#'
#' @param annotated_results A data table of results that have been annotated with the PGRM
#' @param LOUD If TRUE then progress info is printed to the terminal. Default TRUE
#'
#' @return An numeric value of the percent of associations that are powered
#'
#' @export
get_powered_rate = function(annotated_results, include_missing_pheno=TRUE,LOUD = TRUE) {
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


#' Compare two annotated result sets
#'
#' @param results1 A data.table of results that have been annotated with the PGRM
#' @param results2 A data.table of results that have been annotated with the PGRM
#'
#' @return foobar
#'
#' @export
#'
#' @param pheno A data.table of phecode phenotypes. Must have columns `person_id` and `phecode`, as
#'   well as `N` which specifies the number of times the phecode occurred for that person_id
#' @param demos A data.table having one row per person in the cohort. Must have
#'   a column `person_id`. May include other columns used for covariates in later function calls.
#' @param phecode A string specifying the phecode
#' @param MCC The minimum code count required for a case (default MCC=2)
#' @param use_exclude_ranges If TRUE then exclude ranges are applied to controls
#' @param check_sex If TRUE then individuals with sex == F in the demos table will be excluded
#'   from male-specific phenotypes, and vice-versa. Sex specific phecodes are specified in
#'   phecode_info. For this function to work, the demos table must have column `sex` and the values
#'   must be 'M' for Male and 'F' for Female
#'
#' @return A data.table with columns `person_id`, `phecode`, `N` (number of phecode instances),
#'   and `pheno` indicating case, control, or exclude status (represented as 1, 0, NA respectively).
#'   Any additional columns provided in demos will also be returned
#'   prevalence.
#'
#' @export

compare_annotated_results = function(results1, results2){
  summary=data.table()
  results1$dataset='results1'
  results2$dataset='results2'
  r_long=rbind(results1[,c('assoc_ID','odds_ratio','P','L95','U95','rep','powered','Power','rOR','rL95','rU95','dataset')],
               results2[,c('assoc_ID','odds_ratio','P','L95','U95','rep','powered','Power','rOR','rL95','rU95','dataset')])
  r_long=data.table(r_long)
  ## Compare RR all
  #m=glm(rep~dataset,data=r_long,family="binomial")
  #P=summary(m)$coeff['datasetresults2','Pr(>|z|)']
  #OR=exp(summary(m)$coeff['datasetresults2','Estimate'])
  #CIs=exp(confint.default(m))
  #L95=CIs['datasetresults2',1]
  #U95=CIs['datasetresults2',2]
  f=fisher.test(table(r_long$rep,r_long$dataset))
  P=f$p.value
  L95=f$conf.int[1]
  U95=f$conf.int[2]
  OR=f$estimate
  print(glue('\n--------------------------------------'))
  print(glue('Replication rate (all) comparison'))
  RR_ALL_r1=get_RR(results1,include="all",LOUD=FALSE)
  RR_ALL_r2=get_RR(results2,include="all",LOUD=FALSE)
  print(glue('Result1 replication rate (overall) = {r1_RR}',
             r1_RR = sprintf('%1.1f%%', 100 * RR_ALL_r1)))
  print(glue('Result2 replication rate (overall) = {r2_RR}',
             r2_RR = sprintf('%1.1f%%', 100 * RR_ALL_r2)))
  print(glue('Logistic regression rep~dataset'))
  print(glue('Odds ratio (95% CIs) {round(OR,4)} ({round(L95,4)} to {round(U95,4)})'))
  print(glue('P-value {P}'))

  ## Compare Power
  #m=glm(powered~dataset,data=r_long,family="binomial")
  #P=summary(m)$coeff['datasetresults2','Pr(>|z|)']
  #OR=exp(summary(m)$coeff['datasetresults2','Estimate'])
  #CIs=exp(confint.default(m))
  #L95=CIs['datasetresults2',1]
  #U95=CIs['datasetresults2',2]
  f=fisher.test(table(r_long[powered==1]$rep,r_long[powered==1]$dataset))
  P=f$p.value
  L95=f$conf.int[1]
  U95=f$conf.int[2]
  OR=f$estimate
  print(glue('\n--------------------------------------'))
  print(glue('Powered comparison'))
  Power_r1=get_powered_rate(results1,LOUD=FALSE)
  Power_r2=get_powered_rate(results2,LOUD=FALSE)
  print(glue('Result1 replication rate (overall) = {r1_power}',
             r1_power = sprintf('%1.1f%%', 100 * Power_r1)))
  print(glue('Result2 replication rate (overall) = {r2_power}',
             r2_power = sprintf('%1.1f%%', 100 * Power_r2)))
  print(glue('Odds ratio (95% CIs) {round(OR,4)} ({round(L95,4)} to {round(U95,4)})'))
  print(glue('P-value {P}'))

  ## Compare RR powered
  m=glm(rep~dataset,data=r_long[powered==1],family="binomial")
  P=summary(m)$coeff['datasetresults2','Pr(>|z|)']
  OR=exp(summary(m)$coeff['datasetresults2','Estimate'])
  CIs=exp(confint.default(m))
  L95=CIs['datasetresults2',1]
  U95=CIs['datasetresults2',2]
  print(glue('\n--------------------------------------'))
  print(glue('Replication rate (Powered) comparison'))
  RR_ALL_r1=get_RR(results1,LOUD=FALSE)
  RR_ALL_r2=get_RR(results2,LOUD=FALSE)
  print(glue('Result1 replication rate (powered) = {r1_RR}',
             r1_RR = sprintf('%1.1f%%', 100 * RR_ALL_r1)))
  print(glue('Result2 replication rate (powered) = {r2_RR}',
             r2_RR = sprintf('%1.1f%%', 100 * RR_ALL_r2)))
  print(glue('Odds ratio (95% CIs) {round(OR,4)} ({round(L95,4)} to {round(U95,4)})'))
  print(glue('P-value {P}'))

  ## Compare RR, control for power
  m=glm(rep~dataset+Power,data=r_long,family="binomial")
  P=summary(m)$coeff['datasetresults2','Pr(>|z|)']
  OR=exp(summary(m)$coeff['datasetresults2','Estimate'])
  CIs=exp(confint.default(m))
  L95=CIs['datasetresults2',1]
  U95=CIs['datasetresults2',2]
  print(glue('\n--------------------------------------'))
  print(glue('Replication comparison, controlling for Power'))
  AER_r1=get_AER(results1,LOUD=FALSE)
  AER_r2=get_AER(results2,LOUD=FALSE)
  print(glue('Result1 Actual:Expected = {round(AER_r1,3)}'))
  print(glue('Result2 Actual:Expected = {round(AER_r2,3)}'))
  print(glue('Odds ratio (95% CIs) {round(OR,4)} ({round(L95,4)} to {round(U95,4)})'))
  print(glue('P-value {P}'))

  r1=results1[,c('assoc_ID','rsID','phecode','phecode_string','category_string','odds_ratio','P','L95','U95','rep','powered','Power','rOR','rL95','rU95')]
  r2=results2[,c('assoc_ID','odds_ratio','P','L95','U95','rep','powered','Power','rOR','rL95','rU95')]
  r=merge(r1,r2,by="assoc_ID")
  overlapping_rows=nrow(r)
  if(overlapping_rows==0){
    print("no overlapping rows in datasets")
    return(0)
  }
  both_powered=nrow(r[powered.x==1 & powered.y==1])
  both_rep=nrow(r[rep.x==1 & rep.y==1])
  both_powered_rep=nrow(r[rep.x==1 & rep.y==1 & powered.x==1 & powered.y==1])

  r$both_powered = 0
  r[powered.x==1 & powered.y==1]$both_powered = 1

  r$both_rep = 0
  r[rep.x==1 & rep.y==1]$both_rep = 1

  r$both_powered_rep = 0
  r[rep.x==1 & rep.y==1 & powered.x==1 & powered.y==1]$both_powered_rep = 1

  table(r$both_powered,r$both_rep)

  print(glue('\n--------------------------------------'))
  print(glue('Compare odds ratios, all'))

  t=t.test(r$rOR.x,r$rOR.y,paired=T)
  print(glue('Odds ratio comparison (all), n={overlapping_rows}'))
  print(glue('P-value {t$p.value}'))
  print(glue('Mean difference {round(t$estimate,4)}'))
  print(glue('95% CI {round(t$conf.int[1],4)} to {round(t$conf.int[2],4)}'))

  print(glue('\n--------------------------------------'))
  print(glue('Compare odds ratios, Powered'))

  t=t.test(r[rep.x==1 & rep.y==1]$rOR.x,r[rep.x==1 & rep.y==1]$rOR.y,paired=T)
  print(glue('Odds ratio comparison (both replicated), n={both_rep}'))
  print(glue('P-value {t$p.value}'))
  print(glue('Mean difference {round(t$estimate,4)}'))
  print(glue('95% CI {round(t$conf.int[1],4)} to {round(t$conf.int[2],4)}'))}


get_pheno = function(pheno, demos ,phecode,MCC=2,use_exclude_ranges=TRUE,check_sex=FALSE){
  checkPhecodeTable(pheno)
  checkDemosTable(demos)
  checkPhecode(phecode)
  checkMCC(MCC)

  phecode_num = N = sex = NULL

  cur_phecode = phecode

  p=pheno[pheno$phecode==cur_phecode,c('person_id','N')]
  p=data.table(p,key="person_id")
  p$pheno = -9
  p[p$N>=MCC]$pheno = 1

  if(use_exclude_ranges==T){
    to_exclude = c()
    pheno$phecode_num = as.numeric(pheno$phecode)
    ranges = pgrm::exclude_ranges[phecode==cur_phecode]
    for(i in 1:nrow(ranges)){
      to_exclude = union(to_exclude,unique(pheno[phecode_num >= ranges[i,]$range_start & phecode_num < ranges[i,]$range_end]$person_id))
    }
    to_exclude=subset(to_exclude, !to_exclude %in% p$person_id)
    to_exclude=data.table(person_person_id=to_exclude)
    names(to_exclude)[1]="person_id"
    to_exclude$N=0
    to_exclude$pheno=-9
    p=rbind(p,to_exclude)
  }
  p=merge(demos, p,by="person_id",all.x=T)
  p[is.na(pheno)]$pheno=0
  p[p$pheno==-9]$pheno = NA
  p[is.na(N)]$N=0

  if(check_sex==TRUE){
    phecode_sex = sex_check_phecode(cur_phecode)
    if(phecode_sex %in% c("F","M")){
      p[sex!=phecode_sex]$pheno=NA
    }
  }
  return(p)
}

#' Run associations in the pgrm on data from a test cohort. Function requires row level phenotype,
#' genotype, and covariate information.
#'
#' @param geno A matrix or 'BEDMatrix' object containing genetic data
#' @param pheno A data.table having one row per person in the cohort. Must have
#' @param demos A data.table having one row per person in the cohort. Must have
#'   a column `person_id`. May include other columns used for covariates in later function calls.
#' @param covariates A list of covariates to be used in the logistic regression. All covariates
#'   listed must be present in the demos table
#' @param PGRM A data.table of the PGRM, generated with get_PGRM()
#' @param MCC An integer specifying the minimum code count required for a case (default MCC=2)
#' @param minimum_case_count An integer specifiying the minimum number of cases required in the test
#'   cohort to be included in the analysis (default minimum_case_count=100)
#' @param use_exclude_ranges If TRUE then exclude ranges are applied to controls
#' @param check_sex If TRUE then individuals with sex == F in the demos table will be excluded
#'   from male-specific phenotypes, and vice-versa. Sex specific phecodes are specified in
#'   phecode_info. For this function to work, the demos table must have column `sex` and the values
#'   must be 'M' for Male and 'F' for Female
#' @param LOUD If TRUE then progress info is printed to the terminal. Default TRUE
#'
#' @return A data.table with annotated results from association tests
#'
#' @export
run_PGRM_assoc = function(geno, pheno, demos,covariates, PGRM,MCC=2,minimum_case_count=100,use_exclude_ranges=TRUE,check_sex=FALSE,LOUD=TRUE){
  person_id = id = N = . = phecode = cases = NULL
  ## Add check PGRM
  checkGenotypes(geno)
  checkPhecodeTable(pheno)
  checkDemosTable(demos)
  checkCovarList(covariates,demos)
  checkMCC(MCC)
  checkMCC(minimum_case_count)
  checkBool(use_exclude_ranges)
  checkBool(check_sex)
  checkBool(LOUD)
  # ## create formulas for glm using covariates; formula_string_no_sex is for sex-specific phenotypes
  covar_list = paste(covariates, collapse ='+')
  formula_string=paste("pheno~genotype+", covar_list, sep='')
  formula_string_no_sex=gsub("\\+sex","",formula_string)
  formula_string_no_sex=gsub("sex\\+","",formula_string_no_sex)

  ## Change person_id to character to ensure merging works
  demos$person_id=as.character(demos$person_id)
  pheno$person_id=as.character(pheno$person_id)

  # filter demos, pheno, and geno for intersection of IDs
  IDs=union(geno@ped$id,demos$person_id)
  demos =demos[person_id %in% IDs]
  #pheno=pheno[person_id %in% IDs]
  geno=select.inds(geno, id %in% IDs)
  pheno =pheno[person_id %in% IDs]

  ## get available phecodes list
  phecode_counts=pheno[N>=MCC, .(cases = .N), by=phecode]
  available_phecodes = phecode_counts[cases>=minimum_case_count]$phecode

  # filter PGRM for available SNPs and phecodes
  SNPs=geno@snps$id
  PGRM=PGRM[PGRM$SNP %in% SNPs & PGRM$phecode %in% available_phecodes,]

  assoc_to_run=unique(PGRM[,c("SNP", "phecode")])
  assoc_num = nrow(assoc_to_run)
  if(assoc_num == 0){
    print("There are no eligable associations to run")
    return(0)}
  if(LOUD==TRUE){
    print(glue('Attempting {assoc_num} association tests'))
    print(glue('Formula: {formula_string}'))}

  results=data.frame()

  for(i in 1:nrow(assoc_to_run)){
    cur_SNP = assoc_to_run[i,]$SNP
    cur_phecode = assoc_to_run[i,]$phecode
    cur_SNP_index=which(SNPs %in% c(cur_SNP))

    if(LOUD==TRUE){
      print(glue('{i} SNP: {cur_SNP} Phecode: {cur_phecode}'))
    }

    g=data.table(as.matrix(geno[,cur_SNP_index]),keep.rownames=TRUE)
    names(g)[1]='person_id'
    names(g)[2]="genotype"
    setkey(g,"person_id")

    g$genotype = abs(g$genotype-2) ## gaston codes things "backwards" from plink. 2== HOM for ref. Flip this around

    d=merge(g,demos,by="person_id")
    d=get_pheno(pheno=pheno,demos=d, phecode=cur_phecode,MCC=MCC,
                use_exclude_ranges=use_exclude_ranges, check_sex=check_sex)

    formula=as.formula(formula_string)
    if(check_sex==TRUE){
      phecode_sex=sex_check_phecode(cur_phecode)
      if(phecode_sex !="B"){
        formula=
        formula=as.formula(formula_string_no_sex)
      }
    }

    n_case=nrow(d[pheno==1,])
    n_control=nrow(d[pheno==0,])
    if(n_case < minimum_case_count ) {
      next}

    m=glm(formula, data=d,family="binomial")
    conf=confint.default(m)
    P=summary(m)$coeff[2,4]
    odds_ratio=exp(summary(m)$coeff[2,1])

    L95=exp(conf[2,1])
    U95=exp(conf[2,2])
    result = data.frame(SNP=cur_SNP, phecode=cur_phecode,cases=n_case,controls=n_control, P=P, odds_ratio=odds_ratio, L95=L95, U95=U95)
    results = rbind(results, result)}
  return(results)
}

#' Generate a genetic risk score (GRS) based on variants in the PGRM.
#'
#' @param PGRM A data.table of the PGRM, generated with get_PGRM()
#' @param geno A matrix or 'BEDMatrix' object containing genetic data
#' @param phecode A string specifying the phecode
#' @param prune If TRUE, then the SNP list is pruded at the R2 threshold specified
#' @param R2 A numeric value between 0 and 1 indicating the desired R2 value
#' @param LOUD If TRUE then progress info is printed to the terminal. Default TRUE
#'
#' @return A data.table with columns `person_id` and `GRS`
#'
#' @export

make_GRS = function(PGRM,geno,phecode,prune=TRUE,R2=0.2,LOUD=TRUE){
  SNP = id = state = eff = ID = . = NULL
  ## Add check PGRM
  checkGenotypes(geno)
  checkPhecode(phecode)
  ## add check bool
  checkBool(prune)
  ## add R2

  cur_phecode=phecode
  print(glue('Doing {cur_phecode}'))
  GRS=data.frame()

  sub_PGRM=PGRM[phecode==cur_phecode]
  sub_PGRM$Beta = log(sub_PGRM$cat_L95)
  #sub_PGRM$Beta = log(sub_PGRM$cat_OR)

  ## try with confidence interval and test difference
  sub_PGRM[sub_PGRM$risk_allele_dir=="ref",]$Beta = sub_PGRM[sub_PGRM$risk_allele_dir=="ref",]$Beta * -1

  sub_PGRM=sub_PGRM[phecode==cur_phecode & SNP %in% geno@snps$id,c("SNP","assoc_ID","Beta","risk_allele_dir")]
  if(prune){
    total_SNP = nrow(sub_PGRM)
  #   #SNPs_for_GRS=get_SNPs_for_GRS(PGRM,geno,cur_phecode,R2=R2)
    to_prune=get_pruned_SNPs(PGRM, geno, cur_phecode, R2)
    sub_PGRM=sub_PGRM[!SNP %in% to_prune]
    SNPs_left = nrow(sub_PGRM)
    pruned_SNP = total_SNP - SNPs_left
    if(LOUD==TRUE){
      print(glue('Pruned {pruned_SNP} for phecode {cur_phecode} with {SNPs_left} SNPs left'))}
  }

  n_SNPs = nrow(sub_PGRM)
  if(n_SNPs<3){
    print(glue('Not enough eligable SNPs. Only {n_SNPs} for phenotype {cur_phecode}'))
    return(GRS)}
  if(LOUD==TRUE){
    print(glue('There are {n_SNPs} SNPs available for {cur_phecode} GRS'))}

  SNPs_to_include=unique(sub_PGRM$SNP)
  sub_geno=select.snps(geno, id %in% SNPs_to_include)
  sub_geno=data.frame(as.matrix(sub_geno),check.names=F,stringsAsFactors = F)
  sub_geno$ID = row.names(sub_geno)

  n_SNP=ncol(sub_geno)-1
  sub_geno=data.table(sub_geno,key="ID")
  sub_geno=melt(sub_geno,id=ncol(sub_geno),measure=1:n_SNP)

  names(sub_geno)[2] = 'SNP'
  names(sub_geno)[3] = 'state'
  sub_geno$state = abs(sub_geno$state-2) ## gaston codes things "backwards" from plink. 2== HOM for ref. Flip this around
  foo=sub_geno[!is.na(state), .( mean_state = mean(state)), by=SNP]
  sub_geno$state=as.numeric(sub_geno$state)
  sub_geno=merge(sub_geno,foo,by='SNP')
  ## set missing to zero
  sub_geno[is.na(state)]$state = sub_geno[is.na(state)]$mean_state

  sub_geno=data.table(sub_geno,key='SNP')

  GRS=merge(sub_geno,sub_PGRM,by='SNP')

  GRS$eff = GRS$state * GRS$Beta

  GRS=GRS[, .(GRS=sum(eff)), by=ID]
  GRS$GRS_SNP_count = n_SNPs
  names(GRS)[1]='person_id'

  return(GRS)
}
