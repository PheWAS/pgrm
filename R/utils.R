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

checkForCIs = function(results) {
  assertNames(names(results), must.include = c('L95','U95'))
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
