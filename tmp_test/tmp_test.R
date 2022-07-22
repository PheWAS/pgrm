# library(data.table)
# library(gaston)
# library(pgrm)
# library(glue)
# geno_file="~/Dropbox (VUMC)/PGRM/code/data/PGRM_MEGA4_Eur_adult_V1"
# covar_file = "~/Dropbox (VUMC)/PGRM/code/data/covariates_MEGA4_EUR.csv"
# pheno_file = "~/Dropbox (VUMC)/PGRM/code/data/phecodes_V1_2_UP.csv"
#
#
# ## Load PGRM
# PGRM=get_PGRM(ancestry="EUR",build="hg19")
#
# cohort_info = read.csv(file="~/Dropbox (VUMC)/PGRM/code/PGRM/PGRM_cohort.csv",stringsAsFactors = F)
#
# ## Read genotype file
# geno=read.bed.matrix(geno_file)
# geno=select.snps(geno, maf > 0.01)
# geno=select.snps(geno, callrate > 0.98)
# SNPs=geno@snps$id
# length(SNPs) ## 5721
# IDs=geno@ped$id
# length(IDs)  # 65682
#
# ## Read covariate file
# covar=read.csv(covar_file,header=TRUE,stringsAsFactors = F)
# names(covar)[1]="person_id"
# covar=data.table(covar,key="person_id")
#
# ages=read.csv("~/Dropbox (VUMC)/PGRM/code/data/covariates_MEGA4_ages.csv",header=TRUE,stringsAsFactors = F)
# names(ages)[1]="person_id"
# ages=data.table(ages,key="person_id")
# nrow(ages)  ## 93726
# covar=merge(ages,covar,by="person_id")
# covar=covar[person_id %in% IDs]
# names(covar)[3]="last_age"
# covar=covar[last_age>=18*365.25]
#
# ## read phenotype file
# pheno=read.csv(pheno_file,header=TRUE,colClasses=c("character","character","integer"),stringsAsFactors = F)
# names(pheno)[1]="person_id"
# names(pheno)[2]="phecode"
# pheno=pheno[pheno$person_id %in% IDs,] ## include only individuals with genotype data
# pheno=data.table(pheno,key="person_id")
#
# demos=covar
# covariates = c('last_age','sex')
# MCC=2
# minimum_case_count=100
# use_exclude_ranges=TRUE
# LOUD=TRUE
# check_sex=TRUE
#
# r=run_PGRM_assoc(geno=geno, pheno=pheno,demos=demos,covariates = c('last_age','sex'),
#                PGRM=PGRM, MCC=2,minimum_case_count=20000,use_exclude_ranges=TRUE,LOUD=TRUE)
#
# r
# PGRM_ALL[phecode=="401"]
#
#
#
# #run_PGRM_assoc = function(geno, pheno, demos,covariates, PGRM,MCC=2,minimum_case_count=100,use_exclude_ranges=TRUE,check_sex=FALSE,LOUD=TRUE){
#
# sex_check_phecode = function(phecode){
#   if(phecode %in% pgrm::phecode_info[sex=="Male"]$phecode){
#     return("M")
#   }
#   if (phecode %in% pgrm::phecode_info[sex=="Female"]$phecode){
#     return("F")
#   }
#   return("B")}
#
#   # ## create formulas for glm using covariates; formula_string_no_sex is for sex-specific phenotypes
#   covar_list = paste(covariates, collapse ='+')
#   formula_string=paste("pheno~genotype+", covar_list, collapse='')
#   formula_string_no_sex=gsub("sex\\+","",formula_string)
#   if(LOUD==TRUE){
#     print(glue('Formula: {formula_string}'))}
#   formula_string=as.formula(formula_string)
#   formula_string_no_sex=as.formula(formula_string_no_sex)
#
#   ## get available phecodes list
#   phecode_counts=pheno[N>=MCC, .(cases = .N), by=phecode]
#   available_phecodes = phecode_counts[cases>=minimum_case_count]$phecode
#
#   # filter PGRM for available SNPs and phecodes
#   SNPs=geno@snps$id
#   PGRM=PGRM[PGRM$SNP %in% SNPs & PGRM$phecode %in% available_phecodes,]
#
#   assoc_to_run=unique(PGRM[,c("SNP", "phecode")])
#   assoc_num = nrow(assoc_to_run)
#   if(LOUD==TRUE){
#     print(glue('Running {assoc_num} associations'))
#   }
#
#   # filter covar and geno for intersection of IDs
#   IDs=union(geno@ped$id,demos$person_id)
#   demos =demos[person_id %in% IDs]
#   #pheno=pheno[person_id %in% IDs]
#   geno=select.inds(geno, id %in% IDs)
#
#   ## Change person_id to character to ensure merging works
#   demos$person_id=as.character(demos$person_id)
#   pheno$person_id=as.character(pheno$person_id)
#
#   results=data.frame()
#
#   for(i in 1:nrow(assoc_to_run)){
#     cur_SNP = assoc_to_run[i,]$SNP
#     cur_phecode = assoc_to_run[i,]$phecode
#     cur_SNP_index=which(SNPs %in% c(cur_SNP))
#
#     if(LOUD==TRUE){
#       print(glue('{i} SNP: {cur_SNP} Phecode: {cur_phecode}'))
#     }
#
#
#     g=data.table(as.matrix(geno[,cur_SNP_index]),keep.rownames=TRUE)
#     names(g)[1]='person_id'
#     names(g)[2]="genotype"
#     setkey(g,"person_id")
#
#     g$genotype = abs(g$genotype-2) ## gaston codes things "backwards" from plink. 2== HOM for ref. Flip this around
#     demos$person_id=as.character(demos$person_id)
#     d=merge(g,demos,by="person_id")
#     d=get_pheno(pheno=pheno,demos=d, phecode=cur_phecode,MCC=MCC,
#                 use_exclude_ranges=use_exclude_ranges, check_sex=check_sex)
#
#     if(check_sex==TRUE){
#       phecode_sex=sex_check_phecode(cur_phecode)
#       if(phecode_sex !="B"){
#         formula_string=formula_string_no_sex
#       }
#     }
#
#     n_case=nrow(d[pheno==1,])
#     n_control=nrow(d[pheno==0,])
#     if(n_case < minimum_case_count ) {
#       next}
#
#     m=glm(formula_string, data=d,family="binomial")
#     #  conf=confint.default(m)
#     #  P=summary(m)$coeff[2,4]
#     #  odds_ratio=exp(summary(m)$coeff[2,1])
#
#     #  L95=exp(conf[2,1])
#     #  U95=exp(conf[2,2])
#     #  result = data.frame(SNP=cur_SNP, phecode=cur_phecode,cases=n_case,controls=n_control, P=P, odds_ratio=odds_ratio, L95=L95, U95=U95)
#     #  results = rbind(results, result)
#   }
#   return(results)
# }
#
