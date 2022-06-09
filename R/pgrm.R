#' @import checkmate
#' @import data.table
#' @import gaston
#' @import genpwr
#' @import reshape2
NULL

open_PGRM = function(ancestry="all",build="hg19",phecode_version="V1.2"){
  ancestry=toupper(ancestry)
  build=tolower(build)
  PGRM_file = "~/Dropbox (VUMC)/PGRM/code/PGRM/PGRM.csv"
  PGRM_AF_file = "~/Dropbox (VUMC)/PGRM/code/PGRM/PGRM_AF.csv"
  print("Open " %++% PGRM_file %++% " Build " %++% build)
  PGRM=read.table(PGRM_file,sep=",",header = T, stringsAsFactors=F,na.strings=c("-"),colClasses=c('numeric',  'character','character', 'character', 'character', 'character','character', 'character', 'character', 'numeric', 'numeric','numeric', 'numeric', 'character', 'numeric','numeric', 'numeric', 'numeric', 'numeric', 'numeric',  'character'))
  PGRM_AF=read.table(file=PGRM_AF_file,sep=",",header = T, stringsAsFactors=F)

  freq_col_name = ancestry %++% "_freq"
  PGRM=merge(PGRM,PGRM_AF[,c("SNP_hg19",freq_col_name)],by="SNP_hg19")

  names(PGRM)[ncol(PGRM)]="AF"

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
    PGRM$ancestry=NULL
  }

  cases_needed_col_name = 'cases_needed_' %++% ancestry
  PGRM$cases_needed = PGRM[cases_needed_col_name]
  PGRM=data.table(PGRM,key=c("SNP","phecode"))
  PGRM[, c("cases_needed_AFR","cases_needed_EAS","cases_needed_EUR","cases_needed_AMR","cases_needed_SAS","cases_needed_ALL"):=NULL]

  return(PGRM)
}
