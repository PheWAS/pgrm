#' @import checkmate
#' @import data.table
#' @import gaston
#' @import genpwr
#' @importFrom reshape2 melt
#' @importFrom DescTools %c%
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
   ancestry=toupper(ancestry)
   build=tolower(build)
   PGRM=pgrm::PGRM_ALL
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
   setnames(PGRM, freq_col_name, "AF")
   cases_needed_col_name = 'cases_needed_' %c% ancestry
   setnames(PGRM, cases_needed_col_name, "cases_needed")

   cols = c("AFR_freq","EAS_freq","EUR_freq","AMR_freq","SAS_freq","ALL_freq","cases_needed_AFR","cases_needed_EAS","cases_needed_EUR","cases_needed_AMR","cases_needed_SAS","cases_needed_ALL")
   cols=cols[!cols %in% c(freq_col_name,cases_needed_col_name)]
   PGRM[, c(cols):=NULL]

   return(PGRM)
}
