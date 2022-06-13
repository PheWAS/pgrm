checkBuild = function(build) {
  if (!build %in% c('hg19','hg38')) {
    stop('Build must equal hg19 or hg38.')
  }
  invisible()}

checkAncestry = function(ancestry) {
  if (!ancestry %in% c('EAS', 'EUR', 'AFR', 'SAS', 'AMR', 'ALL')) {
    stop('Ancestry must be EAS, EUR, AFR, SAS, AMR, or ALL.')
  }
  invisible()}

checkPhecodeVersion = function(phecode_version) {
  if (!phecode_version %in% c('V1.2','X')) {
    stop('phecode_version must be V1.2 or X')
  }
  invisible()}

check_RR_include = function(include){
  if (!include %in% c('powered','all')) {
    stop('include must be set to powered or all')
  }
  invisible()}

checkResults = function(results) {
  assertDataFrame(results)
  assertNames(names(results), must.include = c('SNP','phecode','cases','controls','odds_ratio','P'))
  invisible()}

checkAnnotatedResults = function(annotated_results,include) {
  assertDataFrame(annotated_results)
  assertNames(names(annotated_results), must.include = c('SNP','phecode','rep'))
  if(include == "powered"){
    assertNames(names(annotated_results), must.include = c('powered'))
  }
  invisible()}

checkAnnotatedResults_forAE = function(annotated_results) {
  assertDataFrame(annotated_results)
  assertNames(names(annotated_results), must.include = c('SNP','phecode','rep','Power'))
  invisible()}
