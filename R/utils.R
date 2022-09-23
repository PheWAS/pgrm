checkBuild = function(build) {
  if (!build %in% c('hg19', 'hg38')) {
    stop('Build must equal hg19 or hg38.')}
  invisible()}


checkAncestry = function(ancestry) {
  if (!ancestry %in% c('EAS', 'EUR', 'AFR', 'SAS', 'AMR', 'ALL')) {
    stop('Ancestry must be EAS, EUR, AFR, SAS, AMR, or ALL')}
  invisible()}


checkPhecodeVersion = function(phecode_version) {
  if (!phecode_version %in% c('V1.2', 'X')) {
    stop('phecode_version must be V1.2 or X')}
  invisible()}



check_RR_include = function(include) {
  if (!include %in% c('powered', 'all')) {
    stop('include must be set to powered or all')}
  invisible()}


checkResults = function(results) {
  assertDataFrame(results)
  assertNames(names(results), must.include = c('SNP', 'phecode', 'cases',
                                               'controls', 'odds_ratio', 'P'))
  invisible()}


checkForCIs = function(results) {
  assertNames(names(results), must.include = c('L95', 'U95'))
  invisible()}

checkBool = function(bool) {
  assertLogical(bool)
  invisible()}

checkSex = function(bool, covars) {
  assertLogical(bool)
  if(bool == TRUE){
    assertDataFrame(covars)
    assertNames(names(covars), must.include = c('person_id', 'sex'))}
  invisible()}

checkR2 = function(R2) {
  assertNumeric(R2, lower = 0, upper = 1)
  invisible()}

checkAnnotatedResults = function(annotated_results, include) {
  assertDataFrame(annotated_results)
  assertNames(names(annotated_results), must.include = c('assoc_ID', 'rep'))
  if (include == 'powered') {
    assertNames(names(annotated_results), must.include = c('powered'))}
  invisible()}

checkIcdTable = function(icds) {
  assertDataFrame(icds)
  assertNames(names(icds), must.include = c('person_id', 'icd','flag'))
  invisible()}


checkAnnotatedResults_forAE = function(annotated_results) {
  assertDataFrame(annotated_results)
  assertNames(names(annotated_results), must.include = c('assoc_ID', 'rep', 'Power'))
  invisible()}

checkPhecodeTable = function(phecode_table) {
  assertDataFrame(phecode_table)
  assertNames(names(phecode_table), must.include = c('person_id', 'phecode', 'N'))
  invisible()}

checkDemosTable = function(demos_table) {
  assertDataFrame(demos_table)
  assertNames(names(demos_table), must.include = c('person_id'))
  invisible()}

checkMCC = function(MCC){
  assertNumeric(MCC, lower=1)}

checkGenotypes = function(genotypes) {
  assertMultiClass(genotypes, c('BEDMatrix', 'matrix','bed.matrix'))
  #assertNames(rownames(genotypes), type = 'unique')
  #assertNames(colnames(genotypes), type = 'unique', disjunct.from = 'score')
  invisible()}

checkCovarList = function(covariate_list, demos_table) {
  assertCharacter(covariate_list)
  demo_cols = names(demos_table)
  if (! all(covariate_list %in% demo_cols)) {
    stop('One or more covariates specified in covariate_list is not in the demos table')
  }
  invisible()}

annotate_power = function(annotated_results, LOUD = FALSE, max_thresh = 50) {
  Power = NULL

  annotated_results[, Power := as.numeric(NA)]
  total = nrow(annotated_results)
  if (isTRUE(LOUD)) {
    print('Doing power calculations')}
  for (i in seq_len(total)) {
    if (isTRUE(LOUD) && i %% 100 == 0) {
      print(glue('{i} of {total}'))}

    annotated_results_tmp = annotated_results[i]

    odds_ratio = annotated_results_tmp$cat_L95
    AF = annotated_results_tmp$RAF
    ## flip AF and OR to minor allele
    if (AF > 0.5) {
      AF = 1 - AF
      odds_ratio = 1 / odds_ratio
    }
    k = annotated_results_tmp$controls / annotated_results_tmp$cases
    N = annotated_results_tmp$controls + annotated_results_tmp$cases
    ## control:case ratio ceiling of max_thresh
    if (k > max_thresh) {
      k = max_thresh
      N = annotated_results_tmp$cases * max_thresh}

    pwr = genpwr.calc(calc = 'power', model = 'logistic', ge.interaction = NULL,
                      Case.Rate = NULL, k = k, N = N, MAF = AF, OR = odds_ratio,
                      Alpha = 0.05, Power = NULL, True.Model = c('Additive'),
                      Test.Model = c('Additive'))
    annotated_results[i]$Power = pwr$Power_at_Alpha_0.05}
  return(annotated_results)}


