checkBuild = function(build) {
  if (!build %in% c('hg19', 'hg38')) {
    stop('Build must equal hg19 or hg38.')}
  invisible()}


checkAncestry = function(ancestry) {
  if (!ancestry %in% c('EAS', 'EUR', 'AFR', 'SAS', 'AMR', 'ALL')) {
    stop('Ancestry must be EAS, EUR, AFR, SAS, AMR, or ALL.')}
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

checkR2 = function(R2) {
  assertNumeric(R2,lower=0,upper=1)
  invisible()}

checkAnnotatedResults = function(annotated_results, include) {
  assertDataFrame(annotated_results)
  assertNames(names(annotated_results), must.include = c('SNP', 'phecode', 'rep'))
  if (include == 'powered') {
    assertNames(names(annotated_results), must.include = c('powered'))}
  invisible()}


checkAnnotatedResults_forAE = function(annotated_results) {
  assertDataFrame(annotated_results)
  assertNames(names(annotated_results), must.include = c('SNP', 'phecode', 'rep', 'Power'))
  invisible()}

checkPhecodeTable = function(phecode_table) {
  assertDataFrame(phecode_table)
  assertNames(names(phecode_table), must.include = c('person_id','phecode','N'))
  invisible()}

checkDemosTable = function(demos_table) {
  assertDataFrame(demos_table)
  assertNames(names(demos_table), must.include = c('person_id'))
  invisible()}

checkMCC = function(MCC){
  assertNumeric(MCC,lower=1)
}

checkPhecode = function(phecode){
  assertCharacter(phecode)
  if (! phecode %in% pgrm::phecode_info$phecode) {
    stop('phecode specified is not a valid phecode')
  }
  #assertTRUE(phecode %in% pgrm::phecode_info$phecode)
}

checkGenotypes = function(genotypes) {
  assertMultiClass(genotypes, c('BEDMatrix', 'matrix','bed.matrix'))
  #assertNames(rownames(genotypes), type = 'unique')
  #assertNames(colnames(genotypes), type = 'unique', disjunct.from = 'score')
  invisible()}

checkCovarList = function(covariate_list,demos_table) {
  assertCharacter(covariate_list)
  demo_cols = names(demos_table)
  if (! all(covariate_list %in% demo_cols)) {
    stop('One or more covariates specified in covariate_list is not in the demos table')
  }
  invisible()}

annotate_power = function(annotated_results, LOUD = FALSE) {
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
    AF = annotated_results_tmp$AF
    ## change AF to risk allele
    if (annotated_results_tmp$risk_allele_dir == 'ref') {
      AF = 1 - AF}
    ## flip AF and OR to minor allele
    if (AF > 0.5) {
      AF = 1 - AF
      odds_ratio = 1 / odds_ratio
    }
    k = annotated_results_tmp$controls / annotated_results_tmp$cases
    N = annotated_results_tmp$controls + annotated_results_tmp$cases
    ## control:case ratio ceiling of 40
    if (k > 40) {
      k = 40
      N = annotated_results_tmp$cases * 40
    }

    pwr = genpwr.calc(calc = 'power', model = 'logistic', ge.interaction = NULL,
                      Case.Rate = NULL, k = k, N = N, MAF = AF, OR = odds_ratio,
                      Alpha = 0.05, Power = NULL, True.Model = c('Additive'),
                      Test.Model = c('Additive'))
    annotated_results[i]$Power = pwr$Power_at_Alpha_0.05
  #  assoc_ID = annotated_results_tmp$assoc_ID
  #  x=pwr$Power_at_Alpha_0.05
  #  print(glue('{assoc_ID} {x} k={k} N={N} AF={AF} OR={odds_ratio}'))

  }
  return(annotated_results)}

sex_check_phecode = function(phecode){
  if(phecode %in% pgrm::phecode_info[sex=='Male']$phecode){
    return('M')}
  if (phecode %in% pgrm::phecode_info[sex=='Female']$phecode){
    return('F')}
  return('B')}

get_pruned_SNPs = function(PGRM, geno, phecode, R2){
  SNP = prune = id = NULL

  to_prune = c()
  cur_phecode = phecode
  SNPs=geno@snps$id

  sub_PGRM=PGRM[phecode==cur_phecode]
  sub_PGRM$prune = 0
  sub_PGRM[!SNP %in% SNPs]$prune = 1 ## mark SNPs that aren't in geno as pruned
  SNP_list=sub_PGRM[prune==0]$SNP

  if(length(SNP_list) <3 ){
    return(to_prune)
  }

  sub_geno=select.snps(geno, id %in% SNP_list)
#  sub_geno=geno[,which(geno@snps$id %in% SNP_list)]
  ld <- LD(sub_geno, c(1,ncol(sub_geno)),measure='r2')

  for(i in 1:length(SNP_list)){
    SNP1 = SNP_list[i]
    chr1 = sub_geno@snps[sub_geno@snps$id==SNP1,]$chr
    if(SNP1 %in% to_prune){
      next
    }
    for(j in 1:length(SNP_list)){
      SNP2 = SNP_list[j]
      chr2 = sub_geno@snps[sub_geno@snps$id==SNP2,]$chr
      if(SNP2 %in% to_prune || SNP1 %in% to_prune || j==i || chr1 != chr2){
        next
      }
      if(ld[SNP1,SNP2]>R2){
        if(sub_PGRM[SNP==SNP1]$cat_LOG10_P > sub_PGRM[SNP==SNP2]$cat_LOG10_P){
          #  print("prune " %c% SNP2 %c% " with R2 " %c% ld[SNP1,SNP2] %c% " to " %c% SNP1)
          to_prune=c(to_prune,SNP2)
        } else {
          # print("prune " %c% SNP1 %c% " with R2 " %c% ld[SNP1,SNP2] %c% " to " %c% SNP2)
          to_prune=c(to_prune,SNP1)
        }
      }
    }
  }
  return(to_prune)
}
