#' @import checkmate
#' @import data.table
#' @import gaston
#' @import genpwr
#' @import DescTools
#' @importFrom reshape2 melt
NULL

#' Get instance of PGRM
#'
#' This function generates a PGRM copy with specified ancestry, build, and phecode version
#'
#' @param ancestry A string that indicates the ancestry of the PGRM. Options EAS, EUR, AFR, SAS, AMR, ALL. Default ALL
#' @param build A string indicating the genome reference build. Options hg19, hg37. Default is hg19.
#' @param phecode_version A string indicating the phecode version. Currently only V1.2 is supported, which is the default
#'
#' @return A data.table of the PGRM
#'
#' @eval example1()
#'
#'
#' @export
get_PGRM = function(ancestry="all",build="hg19",phecode_version="V1.2"){

  ## need functions to check arguments
   ancestry=toupper(ancestry)
   build=tolower(build)
   PGRM=PGRM_ALL
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

   #cols = c("AFR_freq","EAS_freq","EUR_freq","AMR_freq","SAS_freq","ALL_freq","cases_needed_AFR","cases_needed_EAS","cases_needed_EUR","cases_needed_AMR","cases_needed_SAS","cases_needed_ALL")
   #cols=cols[!cols %in% c(freq_col_name,cases_needed_col_name)]
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
#' @eval example2()
#'
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
#' @param build A string indicating the genome reference build used in the results table. Options hg19, hg37. Default is hg19.
#' @param phecode_version A string indicating the phecode version used in the results table. Currently only V1.2 is supported, which is the default
#' @param calculate_power If TRUE then power calculations will be conducted using case and control counts from the results file. Necessary for get_AE(). Default FALSE
#' @param annotate_CI_overlap If TRUE then a column called "annotate_CI_overlap" is added to the table, values:
#'   (**overlap**: 95% CIs of PGRM and test cohort overlap, **test_cohort_greater**: 95% CI of test cohort greater than PGRM, **PGRM_greater**: 95% CI of PGRM greater than test cohort)
#' @param LOUD If TRUE then progress info is printed to the terminal. Default FALSE
#'
#' @return A data.table of the results file annotated with columsn from the PGRM
#'
#' @eval example3()
#'
#'
#' @export
annotate_results = function(results, use_allele_dir=T,ancestry="all",build="hg19",phecode_version="V1.2",calculate_power=FALSE,annotate_CI_overlap=T,LOUD=TRUE){
  PGRM=get_PGRM(ancestry=ancestry,build=build)

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

#' Replicaiton rate (RR) of powered associations
#'
#' This function calculates the replicaiton rate in a test cohort
#'
#' @param results An data.table of results, annotated with the pgrm
#' @param LOUD If TRUE then progress info is printed to the terminal. Default TRUE
#'
#' @return The replication rate of the result set at 80% power
#'
#' @eval example4()
#'
#'
#' @export
get_RR_power = function(results,LOUD=TRUE){
  results=data.table(results)
  powered=nrow(results[powered==1])
  powered_and_rep=nrow(results[powered==1 & rep==1])
  RR=powered_and_rep/powered
  print("Replicated " %c% powered_and_rep %c% " of " %c% powered %c% " for RR=" %c% sprintf("%1.1f%%", 100*RR))
  return(RR)
}
