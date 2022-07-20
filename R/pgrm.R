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
#' @details This function simply returns a copy of PGRM_ALL with specified columns for SNP, AF, and cases_needed
#' according to the arguments specified. The function assigns a column "SNP" as either SNP_hg19 or SNP_hg38, depending
#' on the build argument. It also assigns the allele frequency (AF) and cases_needed columns according to the ancestry.
#' If the ancestry is != "ALL", then the funciton also subsets PGRM_ALL to include only associations that were based in
#' cohorts with the specified ancestry. This funciton is called inside annotate_results(),
#' which can be used to annoate the results of a specific test cohort.
#'
#' @seealso [PGRM_ALL], [annotate_results()]
#'
#' @examples
#' library(pgrm)
#'
#' ## Get a copy of the PGRM for build hg19, East Asian ancestry
#' get_PGRM(build="hg19",ancestry="EAS")
#'
#' @export
get_PGRM = function(ancestry="all",build="hg19",phecode_version="V1.2",unique=T){

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
#' Necessary for get_AE(). Default FALSE
#' @param annotate_CI_overlap If TRUE then a column called "annotate_CI_overlap" is added to the table, values:
#' **overlap**: 95% CIs of PGRM and test cohort overlap
#' **test_cohort_greater**: 95% CI of test cohort greater than PGRM
#' **PGRM_greater**: 95% CI of PGRM greater than test cohort
#' If annotate_CI_overlap is TRUE, then results must include 95% CIs
#' @param LOUD If TRUE then progress info is printed to the terminal. Default FALSE
#'
#' @return A data.table of the results file annotated with columsn from the PGRM
#'
#' @details This function takes a dataframe with summary statistics from a test cohort. For an example of
#' the way to format the results data frame, see one of the results sets included in the package (e.g. results_MGI). (NOTE: If the direction
#' of effect is used to determine if an association is replicated, then the odds ratios of the result set must be oriented to the alternative allele.)
#'
#' The function returns a data.table with the following annotations:
#' \itemize{
#'   \item Phecode informtion, including phecode_string and phecode_category
#'   \item Allele frequencies from GnomAD (column AF), ancestry specified by the ancestry argument
#'   \item The rsID
#'   \item The direction of effect (ref or alt) and risk allele of the original association
#'   \item Summary statistics from the GWAS catalog association, including the -log10(P), odds ratio, and 95% confidence intervals (cat_LOG10_P, cat_OR, cat_L95, cat_U95)
#'   \item The study accession ID from the GWAS catalog
#'   \item A column called powered, which is 1 or 0 indicating whether the test association is powered > 80% (1 if cases >= cases_needed)
#'   \item A column called rep that indicates if the association is replicated (i.e. p<0.05 in the test cohort; if use_allele_dir==TRUE, then the direction of effect from the test cohort must also be consistant with what is reported in the catalog)
#'   \item If annotate_CI_overlap is true, then information about the relationship between the 95% CIs from the catalog and the test set is included in column CI_overlap, and new columns
#' for odds_ratio, L95, and U95 are created (rOR, rL95, rU95) that are oriented to the risk allele. (This option requires that the confidence intervals are reported in the test cohort summary statistics)
#' }
#'
#' @examples
#' library(pgrm)
#' ## annotate the UK Biobank results set. (Filter by cohort_match==0 to exclude associations that were based in whole or in part of UKBB cohort)
#' annotated_results = annotate_results(results_UKBB[cohort_match==0],ancestry="EUR", build="hg38",calculate_power=TRUE)
#'
#' ## Get the replication rate of associations powered at >80%
#' get_RR(annotated_results)
#' ## Get the replication rate of all associations (powered or not)
#' get_RR(annotated_results,include="all")
#'
#' ## Get the actual:expected ratio
#' get_AE(annotated_results)
#'
#' @seealso Example result sets: [results_BBJ], [results_UKBB], [results_BioVU_EUR], [results_BioVU_AFR], [results_MGI]
#' Functions that take annotated result sets [get_RR()],  [get_EA()]
#'
#' @export
annotate_results = function(results, use_allele_dir=T,ancestry="all",build="hg19",phecode_version="V1.2",calculate_power=FALSE,annotate_CI_overlap=T,LOUD=TRUE){
  PGRM=get_PGRM(ancestry=ancestry,build=build,phecode_version=phecode_version)
  checkResults(results)
  if(annotate_CI_overlap==TRUE){
    checkForCIs(results)
  }
  results=data.table(results)
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
#' @examples
#' library(pgrm)
#' ## annotate the UK Biobank results set
#' annotated_results = annotate_results(results_BioVU_AFR,ancestry="AFR", build="hg19",calculate_power=TRUE)
#'
#' ## Get the replication rate of associations powered at >80%
#' get_RR(annotated_results)
#' ## Get the replication rate of all associations (powered or not)
#' get_RR(annotated_results,include="all")
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


#' Run association tests for phenotype/genotype pairs in the PGRM
#'
#' This function takes row level data from a test cohort and runs association tests for all phenotype/genotype pairs included in the PGRM
#'
#' @param genotypes A 'BEDMatrix' object linked to the PLINK bed file containing genetic data. The row names correspond to the person_id's in demos and scores tables. The column names correspond to variant IDs.
#' @param pheno A data.table of phenotypes to be used in the association analysis. Must have columns person_id, phecodes
#' @param demos A data.table of covariates to be used in the association analysis. Must have column person_id and columns for every covariate in covariates list.
#' @param PGRM A copy of the PGRM. Can be generated with function get_PGRM()
#' @param MCC The minimum code count needed for cases. Default == 2.
#' @param minimum_case_count The minimum number of cases required for an association test to be conducted.
#' @param use_exclude_range If TRUE then exclude ranges are applied to controls.
#' @param LOUD If TRUE then progress info is printed to the terminal. Default FALSE
#'
#' @return A data.table of association results for all eligable phenotype/genotype pairs in the PGRM.
#'
#' @details This function takes row level data from a test cohort as well as a copy of the PGRM. It runs an association analysis for each phenotype/genotype pair in the PGRM
#' and returns the results.
#'
#' @examples
#' library(pgrm)
#' ## Get a copy of the PGRM for build hg19, East Asian ancestry
#' PGRM_EAS=get_PGRM(build="hg19",ancestry="EAS")
#'
#' ## Read in a phenotype file
#' pheno=read.csv(pheno_file,header=TRUE,colClasses=c("character","character","integer"),stringsAsFactors = F)
#'
#' ## Read in demographics file
#' demos=read.csv(demos_file,header=TRUE,stringsAsFactors = F)
#'
#' ## Read in genotype file
#' genotypes=read.bed.matrix(geno_file)
#' ## Filter out genotypes with low callrate
#' genotypes=select.snps(geno, callrate > 0.98)
#'
#' covariates=c('sex','last_age','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8')
#' results_MCC2 = run_PGRM_assoc(geno=genotypes, pheno=pheno, demos=demos, covariates=covariates,PGRM=PGRM_EAS,MCC=2,minimum_case_count=100,use_exclude_range=T,LOUD=T)
#'
#' @export
run_PGRM_assoc = function(genotypes, pheno, demos,covariates, PGRM,MCC=2,minimum_case_count=100,LOUD=TRUE,use_exclude_ranges=F){

  ## create formulas for glm using covariates; formula_string_no_sex is for sex-specific phenotypes
  formula_string="pheno~genotype+" %c% paste(covariates, collapse ='+')
  formula_string
  formula_string_no_sex=gsub("sex\\+","",formula_string)
  formula=as.formula(formula_string)
  formula_no_sex=as.formula(formula_string_no_sex)

  ## add a check to make sure covariates are in demos file

  ## get phecode counts
  # phecode_counts=pheno[N>=MCC, .(cases = .N), by=phecode]
  # available_phecodes = phecode_counts[cases>=min_case_count]$phecode
  #
  # ## filter PGRM for available SNPs and phecodes
  #
  # PGRM=PGRM[PGRM$SNP %in% SNPs & PGRM$phecode %in% available_phecodes,]
  #
  # assoc_to_run=unique(PGRM[,c("SNP", "phecode")])
  # assoc_num = nrow(assoc_to_run)
  # print("Running " %c% assoc_num %c% " associations")
  # #assoc_num = 44
  # results=data.frame()
  #
  # ## filter demos and geno for intersection of IDs
  # IDs=geno@ped$id
  # demos =demos[ID %in% IDs]
  # geno=select.inds(geno, id %in% demos$ID)
  # IDs=geno@ped$id
  # pheno=pheno[pheno$ID %in% IDs,]
  #
  # for(i in 1:nrow(assoc_to_run)){
  #
  #   cur_SNP = assoc_to_run[i,]$SNP
  #   cur_phecode = assoc_to_run[i,]$phecode
  #   cur_SNP_index=which(SNPs %in% c(cur_SNP))
  #
  #   g=data.frame(as.matrix(geno[,cur_SNP_index]))
  #   g$ID = row.names(g)
  #   names(g)[1]="genotype"
  #   g=data.table(g,key="ID")
  #
  #   g$genotype = abs(g$genotype-2) ## gaston codes things "backwards" from plink. 2== HOM for ref. Flip this around
  #   d=merge(g,demos,by="ID")
  #
  #   p=get_pheno(pheno, cur_phecode,MCC,use_exclude_ranges = use_exclude_ranges)
  #
  #   d=merge(d,p,by="ID",all.x=T)
  #   d[is.na(pheno)]$pheno=0
  #   d=d[pheno!=-9]
  #
  #   phecode_sex=sex_check_phecode(cur_phecode)
  #   if(phecode_sex !="B"){
  #     # print("sex specific phecode")
  #     d=d[sex == phecode_sex]
  #   }
  #
  #   n_case=nrow(d[pheno==1,])
  #   if(n_case < minimum_case_count ) {
  #     next
  #   }
  #   n_control=nrow(d[pheno==0,])
  #   table(d$genotype)
  #   if(phecode_sex =="B"){
  #     m=glm(formula, data=d,family="binomial")
  #   } else {
  #     m=glm(formula_string_no_sex, data=d,family="binomial")
  #   }
  #
  #   conf=confint.default(m)
  #
  #   P=summary(m)$coeff[2,4]
  #   odds_ratio=exp(summary(m)$coeff[2,1])
  #   if(LOUD==TRUE){
  #     print("[" %c% i %c% "] SNP: " %c% cur_SNP %c% " Phecode: " %c% cur_phecode %c% " P: " %c% P)
  #   }
  #
  #   L95=exp(conf[2,1])
  #   U95=exp(conf[2,2])
  #   result = data.frame(SNP=cur_SNP, phecode=cur_phecode,cases=n_case,controls=n_control, P=P, odds_ratio=odds_ratio, L95=L95, U95=U95)
  #   results = rbind(results, result)
  # }
  # return(results)
}
