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
#'
#' @return A data.table of the PGRM
#'
#' @details This function simply returns a filtered version of PGRM_ALL. The PGRM_ALL data.table includes
#' columns for SNPs specified by build hg19 and hg38 (SNP_hg19 and SNP_hg38). This function assigns a column "SNP" as
#' either SNP_hg19 or SNP_hg38, depending on the build argument. It also assigns the allele frequency (AF) and cases_needed
#' columns according to the ancestry. If the ancestry is != "ALL", then the funciton also subsets PGRM_ALL to include only
#' associations that were based in cohorts with the specified ancestry. This funciton is called inside annotate_results(),
#' which can be used to annoate the results of a specific test cohort.
#'
#' @seealso [PGRM_ALL, annotate_results()]
#'
#' @examples
#' library(pgrm)
#'
#' ## Get a copy of the PGRM for build hg19, East Asian ancestry
#' get_PGRM(build="hg19",ancestry="EAS")
#'
#' @export
get_PGRM = function(ancestry="all",build="hg19",phecode_version="V1.2"){

   ancestry=toupper(ancestry)
   build=tolower(build)
   checkBuild(build)
   checkAncestry(ancestry)
   checkPhecodeVersion(phecode_version)

   PGRM=copy(PGRM_ALL)
   if(build=="hg19"){
     PGRM$SNP_hg38=NULL
     names(PGRM)[1] = "SNP"
   }
   if(build=="hg38"){
     PGRM$SNP_hg19=NULL
     names(PGRM)[2] = "SNP"
   }
  if(ancestry != 'ALL'){
    PGRM = PGRM[PGRM$ancestry==ancestry,]
  }

   freq_col_name = ancestry %c% "_freq"
   cases_needed_col_name = 'cases_needed_' %c% ancestry
   setnames(PGRM, freq_col_name, "AF")
   setnames(PGRM, cases_needed_col_name, "cases_needed")
   PGRM=PGRM[, c('assoc_ID','SNP','ancestry', 'rsID','risk_allele_dir','risk_allele','AF','phecode','phecode_string','category_string','cat_LOG10_P','cat_OR','cat_L95','cat_U95','cases_needed','Study_accession')]


   return(PGRM)
}


#' Add power annotations to a result set that's been merged with PGRM
#'
#' This function adds a power calculation to a result set that has been annotated with PGRM
#'
#' @param annotated_results A data.table of results annotated by PGRM
#' @param LOUD If TRUE then progress info is printed to the terminal. Default FALSE
#'
#' @return A data.table of the results with a column Power added which includes 80% power calculations (alpha=0.05)
#'
#' @export
#'
annotate_power = function(annotated_results,LOUD=FALSE){

  annotated_results$Power = NA
  annotated_results$Power = as.numeric(annotated_results$Power)
  total = nrow(annotated_results)
  if(LOUD==TRUE  ) {
    print("Doing power calculations")
  }
  for(i in 1:nrow(annotated_results)){
    if(LOUD==TRUE & i %% 100 == 0) {
      print(i %c% " of " %c% total)
    }
    odds_ratio = annotated_results[i,]$cat_L95
    AF=annotated_results[i,]$AF
    ## change AF to risk allele
    if(annotated_results[i,]$risk_allele_dir == 'ref'){
      AF=1-AF
    }
    ## flip AF and OR to minor allele
    if(AF>.5){
      AF=1-AF
      odds_ratio = 1/odds_ratio
    }
    k = annotated_results[i,]$controls/annotated_results[i,]$cases
    N = annotated_results[i,]$controls+annotated_results[i,]$cases
    ## control:case ratio ceiling of 20
    if(k>20){
      k=20
      N = annotated_results[i,]$cases * 20
    }
    pwr <- genpwr.calc(calc = "power", model = "logistic", ge.interaction = NULL,
                       Case.Rate=NULL, k=k,N=N,
                       MAF=AF, OR=odds_ratio,Alpha=0.05,Power=NULL,
                       True.Model=c("Additive"),  Test.Model=c( "Additive"))
    annotated_results[i,]$Power = pwr$Power_at_Alpha_0.05
  }
  return(annotated_results)
}



#' Annotate a result set with the PGRM
#'
#' This function annotates a result from a test cohort with information from the PGRM
#'
#' @param results A data frame with results of a test cohort; columns for SNP, phecode, cases, controls, odds_ratio, P (see demo files like results_BBJ for example)
#' @param use_allele_dir If TRUE, direction of effect is used when assessing if an association is replicated
#' @param ancestry A string that specifies ancestry of the PGRM that is then used to annotate the results file. Options EAS, EUR, AFR, SAS, AMR, ALL. Default ALL
#' @param build A string indicating the genome reference build used in the results table. Options hg19, hg38. Default is hg19.
#' @param phecode_version A string indicating the phecode version used in the results table. Currently only V1.2 is supported, which is the default
#' @param calculate_power If TRUE then power calculations will be conducted using case and control counts from the results file. Necessary for get_AE(). Default FALSE
#' @param annotate_CI_overlap If TRUE then a column called "annotate_CI_overlap" is added to the table, values:
#'   (**overlap**: 95% CIs of PGRM and test cohort overlap, **test_cohort_greater**: 95% CI of test cohort greater than PGRM, **PGRM_greater**: 95% CI of PGRM greater than test cohort)
#' @param LOUD If TRUE then progress info is printed to the terminal. Default FALSE
#'
#' @return A data.table of the results file annotated with columsn from the PGRM
#'
#' @export
annotate_results = function(results, use_allele_dir=T,ancestry="all",build="hg19",phecode_version="V1.2",calculate_power=FALSE,annotate_CI_overlap=T,LOUD=TRUE){
  PGRM=get_PGRM(ancestry=ancestry,build=build,phecode_version=phecode_version)
  checkResults(results)

  results$SNP = toupper(results$SNP)
  results=merge(results,PGRM,by=c("SNP","phecode"))
  results$powered = 0
  results[!is.na(cases_needed) & cases>=cases_needed,]$powered=1

  results$rep = 0
  results[results$P<0.05,]$rep = 1

  if(use_allele_dir){
    results[risk_allele_dir=="ref" & odds_ratio > 1,]$rep = 0
    results[risk_allele_dir=="alt" & odds_ratio <1,]$rep = 0
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
#' This function calculates the replicaiton rate in a test cohort that has been annotated with PGRM. By default, it calculates the RR for associations that powered at >80%.
#'
#' @param annotated_results A data table of results that have been annotated with the PGRM
#' @param include A character string. If "powered" then only powered associations are included (default). If "all" then all associations are included
#' @param LOUD If TRUE then progress info is printed to the terminal. Default TRUE
#'
#' @return An numeric value of the replication rate of the result set
#'
#' @export
get_RR = function(annotated_results,include="powered",LOUD=TRUE){

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
get_AE = function(annotated_results,LOUD=TRUE){
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
