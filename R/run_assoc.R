#'
#'
#' #' Make a phecode phenotype file from table of ICD codes
#' #'
#' #' This function calculates % of powered associations for a test cohort that has been annotated with PGRM.
#' #'
#' #' @param icds A data table ICD codes. Must include columns `person_id`, `icd`, `flag`
#' #' @param phecode_version Phecode map version. Currently only supports 'V1.2', which is default
#' #'
#' #' @return An numeric value of the percent of associations that are powered
#' #'
#' #' @export
#' make_pheno = function(icds, phecode_version = 'V1.2') {
#'   checkIcdTable(icds)
#'   checkPhecodeVersion(phecode_version)
#'   if(phecode_version == 'V1.2') {
#'     icdPhecodeMap = icdPhecodeMap_V1_2}
#'   pheno = merge(icds, icdPhecodeMap, by = c('icd', 'flag'), allow.cartesian = TRUE)
#'   pheno = pheno[, c('person_id', 'phecode', 'entry_date')]
#'   pheno = unique(pheno)
#'   pheno = pheno[, .N, by = c('person_id', 'phecode')]
#'   return(pheno)}
#'
#' #' Add case/control status for a specified phenotype to a covariate file.
#' #'
#' #' @param pheno A data.table with phecode phenotypes. columns `person_id`, `phecode`, and `N`
#' #' required. The phecode columns should be a character vector. Column `N` specifies the number of
#' #' unique dates a phecode occurs for that individual. Only phecodes that are present for an individual
#' #' are included (i.e. N is always greater than one).
#' #' @param covars A data.table having one row per person in the cohort. Must havea column names
#' #' `person_id`. Must include other columns specified in covariate_names.
#' #' @param phecode A character vector of the phecode
#' #' @param MCC A numeric value specified the minimum code count required for cases.
#' #' @param use_exclude_ranges If TRUE then individuals who are not cases but have a phecode within
#' #' the exclude ranges are excluded.
#' #' @param check_sex If TRUE then the `pheno` column is set to zero for individuals with a sex that
#' #' does not match the specified sex. This argument is only relevant for sex specific codes
#' #' (e.g. Prostate cancer). To use this function, sex must be included as a column in the covars data.table.
#' #'
#' #' @return The covar data.table with a new column called `pheno`. This column is 0 if the individual,
#' #' is a control for specified phecode, 1 if they are a case, and NA if they are excluded.
#' #'
#' #' @export
#' get_pheno = function(pheno, covars, phecode, MCC = 2, use_exclude_ranges = TRUE, check_sex = FALSE) {
#'   checkPhecodeTable(pheno)
#'   checkDemosTable(covars)
#'   checkPhecode(phecode)
#'   checkSex(check_sex, covars)
#'   checkMCC(MCC)
#'
#'   phecode_num = N = sex = NULL
#'
#'   cur_phecode = phecode
#'
#'   p = pheno[pheno$phecode == cur_phecode, c('person_id', 'N')]
#'   p = data.table(p, key = 'person_id')
#'   p$pheno = -9
#'   p[p$N >= MCC]$pheno = 1
#'
#'   if(use_exclude_ranges == TRUE){
#'     to_exclude = c()
#'     pheno$phecode_num = as.numeric(pheno$phecode)
#'     ranges = pgrm::exclude_ranges[phecode == cur_phecode]
#'     for(i in 1:nrow(ranges)){
#'       to_exclude = union(to_exclude, unique(pheno[phecode_num >= ranges[i,]$range_start &
#'                                                     phecode_num < ranges[i,]$range_end]$person_id))
#'     }
#'     to_exclude = subset(to_exclude, !to_exclude %in% p$person_id)
#'     to_exclude = data.table(person_person_id = to_exclude)
#'     names(to_exclude)[1] = 'person_id'
#'     to_exclude$N = 0
#'     to_exclude$pheno = -9
#'     p = rbind(p, to_exclude)
#'   }
#'   p = merge(covars, p, by = 'person_id', all.x=TRUE)
#'   p[is.na(pheno)]$pheno = 0
#'   p[p$pheno == -9]$pheno = NA
#'   p[is.na(N)]$N = 0
#'
#'   if(check_sex == TRUE){
#'     phecode_sex = sex_check_phecode(cur_phecode)
#'     if(phecode_sex %in% c('F', 'M')){
#'       p[sex != phecode_sex]$pheno = NA
#'     }
#'   }
#'   return(p)}
#'
#' #' Run associations in the pgrm on data from a test cohort. Function requires row level phenotype,
#' #' genotype, and covariate information.
#' #'
#' #' (NOTE: If you used principle components as covariates, make sure they are not coded as numeric.
#' #' Otherwise the logistic regression will treat them as categorical variables, making the function
#' #' run very slowly before failing.)
#' #'
#' #' @param geno A matrix or 'BEDMatrix' object containing genetic data
#' #' @param pheno A data.table with phecode phenotypes. columns `person_id`, `phecode`, and `N`
#' #' required. The phecode columns should be a character vector. Column `N` specifies the number of
#' #' unique dates a phecode occurs for that individual. Only phecodes that are present for an individual
#' #' are included (i.e. N is always greater than one).
#' #' @param covars A data.table having one row per person in the cohort. Must havea column names
#' #' `person_id`. Must include other columns specified in covariate_names.
#' #' @param covariate_names A list of covariates to be used in the logistic regression. All covariates
#' #'   listed must be present in the covars table
#' #' @param PGRM A data.table of the PGRM, generated with get_PGRM()
#' #' @param MCC An integer specifying the minimum code count required for a case (default MCC=2)
#' #' @param minimum_case_count An integer specifiying the minimum number of cases required in the test
#' #'   cohort to be included in the analysis (default minimum_case_count=100)
#' #' @param use_exclude_ranges If TRUE then exclude ranges are applied to controls
#' #' @param check_sex If TRUE then individuals with sex == F in the covars table will be excluded
#' #'   from male-specific phenotypes, and vice-versa. Sex specific phecodes are specified in
#' #'   phecode_info. For this function to work, the covars table must have column `sex` and the values
#' #'   must be 'M' for Male and 'F' for Female
#' #' @param LOUD If TRUE then progress info is printed to the terminal. Default TRUE
#' #'
#' #' @return A data.table with annotated results from association tests
#' #'
#' #' @eval run_PGRM_assoc_ex()
#' #'
#' #' @export
#' run_PGRM_assoc = function(geno, pheno, covars,covariate_names, PGRM, MCC=2, minimum_case_count = 100,
#'                           use_exclude_ranges = TRUE, check_sex = FALSE, LOUD = TRUE) {
#'   person_id = id = N = . = phecode = cases = NULL
#'   ## Add check PGRM
#'   checkGenotypes(geno)
#'   checkPhecodeTable(pheno)
#'   checkDemosTable(covars)
#'   checkCovarList(covariate_names,covars)
#'   checkMCC(MCC)
#'   checkMCC(minimum_case_count)
#'   checkBool(use_exclude_ranges)
#'   checkSex(check_sex, covars)
#'   checkBool(LOUD)
#'   # ## create formulas for glm using covariates; formula_string_no_sex is for sex-specific phenotypes
#'   covar_list = paste(covariate_names, collapse = '+')
#'   formula_string = paste('pheno~genotype+', covar_list, sep = '')
#'   formula_string_no_sex = gsub('\\+sex', '', formula_string)
#'   formula_string_no_sex = gsub('sex\\+', '', formula_string_no_sex)
#'
#'   ## Change person_id to character to ensure merging works
#'   covars$person_id = as.character(covars$person_id)
#'   pheno$person_id = as.character(pheno$person_id)
#'
#'   # filter covars, pheno, and geno for intersection of IDs
#'   IDs = union(geno@ped$id, covars$person_id)
#'   covars = covars[person_id %in% IDs]
#'   geno = select.inds(geno, id %in% IDs)
#'   pheno = pheno[person_id %in% IDs]
#'
#'   ## get available phecodes list
#'   phecode_counts = pheno[N >= MCC, .(cases = .N), by = phecode]
#'   available_phecodes = phecode_counts[cases >= minimum_case_count]$phecode
#'
#'   # filter PGRM for available SNPs and phecodes
#'   SNPs = geno@snps$id
#'   PGRM = PGRM[PGRM$SNP %in% SNPs & PGRM$phecode %in% available_phecodes,]
#'
#'   assoc_to_run = unique(PGRM[, c('SNP', 'phecode')])
#'   assoc_num = nrow(assoc_to_run)
#'   if(assoc_num == 0){
#'     print('There are no eligable associations to run')
#'     return(0)}
#'   if(LOUD == TRUE){
#'     print(glue('Attempting {assoc_num} association tests'))
#'     print(glue('Formula: {formula_string}'))}
#'
#'   results = data.frame()
#'
#'   for(i in 1:nrow(assoc_to_run)){
#'     cur_SNP = assoc_to_run[i,]$SNP
#'     cur_phecode = assoc_to_run[i,]$phecode
#'     cur_SNP_index = which(SNPs %in% c(cur_SNP))
#'
#'     if(LOUD == TRUE){
#'       print(glue('{i} SNP: {cur_SNP} Phecode: {cur_phecode}'))
#'     }
#'
#'     g = data.table(as.matrix(geno[, cur_SNP_index]), keep.rownames = TRUE)
#'     names(g)[1] = 'person_id'
#'     names(g)[2] = 'genotype'
#'     setkey(g, 'person_id')
#'
#'     g$genotype = abs(g$genotype - 2) ## gaston codes 2==HOM for ref. Flip this around
#'
#'     d = merge(g, covars, by='person_id')
#'     d = get_pheno(pheno = pheno, covars = d, phecode = cur_phecode, MCC = MCC,
#'                   use_exclude_ranges = use_exclude_ranges, check_sex = check_sex)
#'summary_statistics
#'     formula = as.formula(formula_string)
#'     if(check_sex == TRUE) {
#'       phecode_sex = sex_check_phecode(cur_phecode)
#'       if(phecode_sex != 'B') {
#'         formula = as.formula(formula_string_no_sex)}}
#'
#'     n_case = nrow(d[pheno == 1,])
#'     n_control = nrow(d[pheno == 0,])
#'     if(n_case < minimum_case_count) {
#'       next}
#'
#'     m = glm(formula, data = d, family = 'binomial')
#'     conf = confint.default(m)
#'     P = summary(m)$coeff[2, 4]
#'     odds_ratio = exp(summary(m)$coeff[2, 1])
#'
#'     L95 = exp(conf[2, 1])
#'     U95 = exp(conf[2, 2])
#'     result = data.frame(SNP = cur_SNP, phecode = cur_phecode, cases = n_case, controls = n_control,
#'                         P = P, odds_ratio = odds_ratio, L95 = L95, U95 = U95)
#'     results = rbind(results, result)}
#'   return(results)}
