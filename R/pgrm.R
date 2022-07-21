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
get_pheno = function(pheno, demos ,phecode,MCC=2,use_exclude_ranges=TRUE,check_sex=FALSE){
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
  p[p$pheno==-9]$pheno = NA
  p[is.na(N)]$N=0

  if(check_sex==TRUE){
    phecode_sex = sex_check_phecode(cur_phecode)
    if(phecode_sex %in% c("F","M")){
      p[sex!=phecode]$pheno=NA
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
#' @param LOUD If TRUE then progress info is printed to the terminal. Default TRUE
#'
#' @return A data.table with annotated results from association tests
#'
#' @export
run_PGRM_assoc = function(geno, pheno, demos,covariates, PGRM,MCC=2,minimum_case_count=100,use_exclude_ranges=TRUE,LOUD=TRUE){

  checkGenotypes(geno)
  checkPhecodeTable(pheno)
  checkDemosTable(demos)
  checkCovarList(covariates,demos)
  checkMCC(MCC)
  checkMCC(minimum_case_count)
  # ## create formulas for glm using covariates; formula_string_no_sex is for sex-specific phenotypes
  # formula_string="pheno~genotype+" %c% paste(covariates, collapse ='+')
  # formula_string
  # formula_string_no_sex=gsub("sex\\+","",formula_string)
  # formula=as.formula(formula_string)
  # formula_no_sex=as.formula(formula_string_no_sex)
  #
  # ## get phecode counts
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
  # ## filter covar and geno for intersection of IDs
  # IDs=geno@ped$id
  # covar =covar[ID %in% IDs]
  # geno=select.inds(geno, id %in% covar$ID)
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
  #   d=merge(g,covar,by="ID")
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

