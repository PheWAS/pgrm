#' @import checkmate
#' @import data.table
#' @import gaston
#' @import genpwr
#' @importFrom DescTools %c%
NULL

#' Get instance of PGRM
#'
#' This function generates a PGRM copy with specified ancestry, build, and phecode version
#'
#' @param ancestry A string that indicates the ancestry of the PGRM. Options EAS, EUR, AFR, SAS, AMR, ALL. Default ALL
#' @param build A string indicating the genome reference build. Options hg19, hg38. Default is hg19.
#' @param phecode_version A string indicating the phecode version. Currently only V1.2 is supported, which is the default
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
#' @seealso [PGRM_ALL], [annotate_results()]
#'
#' @eval example_get_PGRM()
#'
#' @export
get_PGRM = function(ancestry="all",build="hg19",phecode_version="V1.2",unique=T){

  ## Avoid warnings about global vars
   cat_LOG10_P = SNP = . = phecode = NULL

   ancestry=toupper(ancestry)
   build=tolower(build)
   checkBuild(build)
   checkAncestry(ancestry)
   checkPhecodeVersion(phecode_version)



   PGRM=copy(pgrm::PGRM_ALL)
   if(build=="hg19"){
     PGRM$SNP_hg38=NULL
     names(PGRM)[1] = "SNP"
   }
   if(build=="hg38"){
     PGRM$SNP_hg19=NULL
     names(PGRM)[2] = "SNP"
   }
  if(ancestry != 'ALL'){
    a = ancestry
    PGRM <- PGRM[ancestry==a]
  }
  if(ancestry == 'ALL' & unique == TRUE) {
    ## make the "ALL" PGRM unique by SNP/phecode
    uniq_PGRM = PGRM[, .(cat_LOG10_P=max(cat_LOG10_P)), by=list(SNP,phecode)]
    PGRM=merge(uniq_PGRM, PGRM,by=c("SNP","phecode","cat_LOG10_P"))
  }
   freq_col_name = ancestry %c% "_freq"
   cases_needed_col_name = 'cases_needed_' %c% ancestry
   setnames(PGRM, freq_col_name, "AF")
   setnames(PGRM, cases_needed_col_name, "cases_needed")
   PGRM <- PGRM[, c('assoc_ID','SNP','ancestry', 'rsID','risk_allele_dir','risk_allele','AF','phecode','phecode_string','category_string','cat_LOG10_P','cat_OR','cat_L95','cat_U95','cases_needed','Study_accession')]


   return(PGRM)
}


#' Annotate a result set with the PGRM
#'
#' This function annotates a result from a test cohort with information from the PGRM
#'
#' @param results A data frame (or data.table) with results of a test cohort; columns for SNP, phecode, cases, controls, odds_ratio,
#' P (see demo files for example (e.g results_MGI))
#' @param use_allele_dir If TRUE, direction of effect is used when assessing if an association is replicated. To use this argument,
#' odds ratios must be reported for the alternative allele
#' @param ancestry A string that specifies ancestry of the PGRM that is then used to annotate the results file.
#' Options EAS, EUR, AFR, SAS, AMR, ALL. Default ALL
#' @param build A string indicating the genome reference build used in the results table. Options hg19, hg38. Default is hg19.
#' @param phecode_version A string indicating the phecode version used in the results table. (Currently only V1.2 is supported)
#' @param calculate_power If TRUE then power calculations will be conducted using case and control counts from the results file.
#' Necessary for get_AER(). Default FALSE
#' @param annotate_CI_overlap If TRUE then a column called "annotate_CI_overlap" is added to the
#' table, values:
#' **overlap**: 95% CIs of PGRM and test cohort overlap
#' **test_cohort_greater**: 95% CI of test cohort greater than PGRM
#' **PGRM_greater**: 95% CI of PGRM greater than test cohort
#' If annotate_CI_overlap is TRUE, then results must include 95% CIs
#' @param LOUD If TRUE then progress info is printed to the terminal. Default FALSE
#'
#' @return A data.table of the results file annotated with columsn from the PGRM
#'
#' @details This function takes a dataframe with summary statistics from a test cohort. For an
#' example of the way to format the results data frame, see one of the results sets included in the
#' package (e.g. results_MGI). (NOTE: If the direction
#' of effect is used to determine if an association is replicated, then the odds ratios of the
#' result set must be oriented to the alternative allele.)
#'
#' The function returns a data.table with the following annotations:
#' \itemize{
#'   \item Phecode informtion, including phecode_string and phecode_category
#'   \item Allele frequencies from GnomAD (column AF), ancestry specified by the ancestry argument
#'   \item The rsID
#'   \item The direction of effect (ref or alt) and risk allele of the original association
#'   \item Summary statistics from the GWAS catalog association, including the -log10(P), odds ratio, and 95% confidence intervals
#'   (cat_LOG10_P, cat_OR, cat_L95, cat_U95)
#'   \item The study accession ID from the GWAS catalog
#'   \item A column called powered, which is 1 or 0 indicating whether the test association is powered > 80% (1 if cases >= cases_needed)
#'   \item A column called rep that indicates if the association is replicated (i.e. p<0.05 in the test cohort; if use_allele_dir==TRUE,
#'   then the direction of effect from the test cohort must also be consistant with what is reported in the catalog)
#'   \item If annotate_CI_overlap is true, then information about the relationship between the 95% CIs from the catalog and the test set is
#'   included in column CI_overlap, and new columns for odds_ratio, L95, and U95 are created
#'   (rOR, rL95, rU95) that are oriented to the risk allele. (This option requires that the confidence
#'   intervals are reported in the test cohort summary statistics)
#' }
#'
#' @eval annotate_UKBB()
#'
#' @seealso Example result sets: [results_BBJ], [results_UKBB], [results_BioVU_EUR], [results_BioVU_AFR], [results_MGI]
#'
#' @export
annotate_results = function(results, use_allele_dir=T,ancestry="all",build="hg19",phecode_version="V1.2",calculate_power=FALSE,annotate_CI_overlap=T,LOUD=TRUE){

  ## Avoid warnings about global vars
  cases = cases_needed = risk_allele_dir = odds_ratio = cat_L95 = cat_U95 = rL95 = rU95 =  NULL

  PGRM=get_PGRM(ancestry=ancestry,build=build,phecode_version=phecode_version)
  checkResults(results)
  if(annotate_CI_overlap==TRUE){
    checkForCIs(results)
  }

  results=data.table(results)
  results$SNP = toupper(results$SNP)
  results=merge(results,PGRM,by=c("SNP","phecode"))
  results$powered = 0



  results[!is.na(cases_needed) & cases>=cases_needed]$powered=1

  results$rep = 0
  results[results$P<0.05,]$rep = 1

  if(use_allele_dir){
    results[risk_allele_dir=='ref' & odds_ratio > 1]$rep = 0
    results[risk_allele_dir=='alt' & odds_ratio <1]$rep = 0
  }
  if(calculate_power==TRUE){
    results=annotate_power(results, LOUD=LOUD)
  }
  if(annotate_CI_overlap==TRUE){
    results$rOR=results$odds_ratio
    results$rL95=results$L95
    results$rU95=results$U95
    results[risk_allele_dir=="ref"]$rOR=1/results[risk_allele_dir=="ref"]$odds_ratio
    results[risk_allele_dir=="ref"]$rU95=1/results[risk_allele_dir=="ref"]$L95
    results[risk_allele_dir=="ref"]$rL95=1/results[risk_allele_dir=="ref"]$U95

    results$CI_overlap = ''
    results[rL95 >= cat_L95 & rL95 <= cat_U95  ]$CI_overlap = 'overlap'
    results[rU95 >= cat_L95  & rL95 <= cat_U95 ]$CI_overlap = 'overlap'
    results[rL95  > cat_U95 ]$CI_overlap = 'test_cohort_greater'
    results[cat_L95  > rU95   ]$CI_overlap = 'PGRM_greater'
  }
  return(results)
}

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
get_RR = function(annotated_results,include="powered",LOUD=TRUE){
  ## Avoid warnings about global vars
  powered = NULL

  include = tolower(include)
  check_RR_include(include)
  checkAnnotatedResults(annotated_results,include)

  if(include == "powered"){
    denominator=nrow(annotated_results[powered==1])
    numerator=nrow(annotated_results[powered==1 & rep==1])
  } else {
    denominator=nrow(annotated_results)
    numerator=nrow(annotated_results[rep==1])
  }
  RR=numerator/denominator
  if(LOUD==TRUE){
    print("Replicated " %c% numerator %c% " of " %c% denominator %c% " for RR=" %c% sprintf("%1.1f%%", 100*RR))
  }
  return(RR)
}

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
get_AER = function(annotated_results,LOUD=TRUE){
  ## Avoid warnings about global vars
  Power = NULL

  checkAnnotatedResults_forAE(annotated_results)
  r=annotated_results[!is.na(Power)]
  expected=sum(r$Power)
  actual=sum(r$rep)
  AE=actual/expected
  expected=round(expected,1)
  total_assoc = nrow(r)
  uniq_phecode=length(unique(r$phecode))
  if(LOUD==TRUE){
  print("Expected " %c% expected %c% ", replicated " %c% actual %c% " for AE=" %c% round(AE,3) %c% " (" %c% total_assoc %c% " associations for " %c% uniq_phecode %c% " uniq phecodes)" )
  }
  return(AE)
}

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
get_powered_rate = function(annotated_results,LOUD=TRUE){
  checkAnnotatedResults(annotated_results, include="powered")
  total_rows=nrow(annotated_results)
  powered=nrow(annotated_results[powered==1])
  powered_rate=powered/total_rows
  pr=round(powered_rate,2)

  if(LOUD==TRUE) {
    uniq_phecode=length(unique(annotated_results[powered==1]$phecode))
    print("Powered for " %c% powered %c% " of " %c% total_rows %c% " associations " %c% pr %c% " for " %c% uniq_phecode %c% " uniq phecodes")
  }
  return(powered_rate)
}

#' Define a phenotype from an ICD file and list of person_id
#'
#' @param pheno A data.table of phecode phenotypes. Must have columns `person_id` and `phecode`, as
#'   well as `N` which specifies the number of times the phecode occurred for that person_id
#' @param demos A data.table having one row per person in the cohort. Must have
#'   a column `person_id`. May include other columns used for covariates in later function calls.
#' @param phecode A string specifying the phecode
#' @param MCC The minimum code count required for a case (default MCC=2)
#' @param use_exclude_ranges If TRUE then exclude ranges are applied to controls
#'
#' @export
get_pheno = function(pheno, demos ,phecode,MCC=2,use_exclude_ranges=FALSE){
  checkPhecodeTable(pheno)
  checkDemosTable(demos)
  checkPhecode(phecode)
  checkMCC(MCC)

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
  p=p[p$pheno!=-9]
  p[is.na(N)]$N=0
  return(p)
}

