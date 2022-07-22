library('data.table')
#registerDoSEQ()

snapshot = function(xObs, path) {
  if (file.exists(path)) {
    xExp = qs::qread(path)
  } else {
    qs::qsave(xObs, path)
    xExp = xObs}
  return(xExp)}

dataDir = 'data'

phecode_table_test = data.table(
  person_id = c(1, rep(2L, 3), 3, rep(4L, 3)),
  phecode = c('250.2', '250.2', '174.11', '153', '250', '250.1', '427.21', '185'),
  N= c(2,2,8,7,8,5,3,7)
)


demos_table_test = data.table(
  person_id = 1:4, sex = c('F', 'M', 'F', 'M'),
  last_age = c(28772,18028,11636,14589))

#phecode_table_test
#get_pheno(pheno=phecode_table_test, demos=demos_table_test, phecode="185",MCC=1,use_exclude_ranges = T,check_sex=T)
#get_pheno(pheno=phecode_table_test, demos=demos_table_test, phecode="185",MCC=1,use_exclude_ranges = T,check_sex=F)

#library(gaston)
#geno=read.bed.matrix('~/Documents/GitHub/pgrm/tests/testthat/data/geno_test')
#geno@snps$id[1]="10:112678657:G:T"
#cur_SNP_index=1
#g=data.table(as.matrix(geno[,cur_SNP_index]),keep.rownames=TRUE)

#run_PGRM_assoc(geno=geno, pheno=phecode_table_test,demos=demos_table_test,covariates = c('last_age','sex'),
#               PGRM=cur_PGRM, MCC=2,minimum_case_count=1,use_exclude_ranges=TRUE,LOUD=TRUE)

#g=data.table(as.matrix(geno[,1]))
#g$person_id = row.names(g)

#head(g)
#names(g)[1]="genotype"
#setkey(g,"person_id")
#head(g)


