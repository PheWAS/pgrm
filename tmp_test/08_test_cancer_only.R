library(data.table)
library(gaston)
library(pgrm)
library(glue)
geno_file="~/Dropbox (VUMC)/PGRM/code/data/PGRM_MEGA4_Eur_adult_V1"
covar_file = "~/Dropbox (VUMC)/PGRM/code/data/covariates_MEGA4_EUR.csv"
pheno_file = "~/Dropbox (VUMC)/PGRM/code/data/phecodes_V1_2_UP.csv"
cancer_reg_file = "~/Dropbox (VUMC)/PGRM/code/data/MEGA_in_cancer_reg.csv"

set.seed(4)

## Load PGRM
PGRM=get_PGRM(ancestry="EUR",build="hg19")

cohort_info = read.csv(file="~/Dropbox (VUMC)/PGRM/code/PGRM/PGRM_cohort.csv",stringsAsFactors = F)
cohort_info=data.table(cohort_info)

PGRM=PGRM[!assoc_ID %in% cohort_info[BioVU_only==1]$assoc_ID]


## Read genotype file
geno=read.bed.matrix(geno_file)
geno=select.snps(geno, maf > 0.01)
geno=select.snps(geno, callrate > 0.98)
SNPs=geno@snps$id
length(SNPs) ## 5721
IDs=geno@ped$id
length(IDs)  # 65682


## Read covariate file
covar=read.csv(covar_file,header=TRUE,stringsAsFactors = F)
names(covar)[1]="person_id"
covar=data.table(covar,key="person_id")

ages=read.csv("~/Dropbox (VUMC)/PGRM/code/data/covariates_MEGA4_ages.csv",header=TRUE,stringsAsFactors = F)
names(ages)[1]="person_id"
ages=data.table(ages,key="person_id")
nrow(ages)  ## 93726
covar=merge(ages,covar,by="person_id")
covar=covar[person_id %in% IDs]
names(covar)[3]="last_age"
covar=covar[last_age>=18*365.25]
covar=covar[UNIQ_DATE>=4]
nrow(covar) ## 62777

## read phenotype file
pheno=read.csv(pheno_file,header=TRUE,colClasses=c("character","character","integer"),stringsAsFactors = F)
names(pheno)[1]="person_id"
names(pheno)[2]="phecode"
pheno=pheno[pheno$person_id %in% IDs,] ## include only individuals with genotype data
pheno=data.table(pheno,key="person_id")

## read cohort info
cancer_reg = read.csv(cancer_reg_file,header=TRUE,stringsAsFactors = F)
nrow(cancer_reg) ## 17696


covar_cancer_reg = covar[covar$person_id %in% cancer_reg$person_id]
pheno_cancer_reg = pheno[pheno$person_id %in% cancer_reg$person_id]
geno_cancer_reg = select.inds(geno, id %in% covar_cancer_reg$person_id)
IDs_cancer_reg=geno_cancer_reg@ped$id
length(IDs_cancer_reg)  # 15302
nrow(covar_cancer_reg) ## 15302



r_cancer_reg=run_PGRM_assoc(geno=geno_cancer_reg, pheno=pheno_cancer_reg,demos=covar_cancer_reg,covariates = c('last_age','sex','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'),
                            PGRM=PGRM, MCC=2,minimum_case_count=100,use_exclude_ranges=TRUE,LOUD=TRUE,check_sex = TRUE)
r_cancer_reg=annotate_results(r_cancer_reg,ancestry="EUR",build="hg19",calculate_power = TRUE)
get_RR(r_cancer_reg) ## Replicated 140 of 208 for RR=67.3%
get_AER(r_cancer_reg) ## Expected 737.3, replicated 544 for AE=0.738 (2785 associations for 74 uniq phecodes)

write.table(r_cancer_reg,file="anno_BioVU_EUR_cancer_reg.csv",row.names = F, col.names = T)



IDs_no_cancer = sample(covar[!covar$person_id %in% cancer_reg$person_id]$person_id,size=nrow(covar_cancer_reg))
covar_no_cancer_reg = covar[covar$person_id %in% IDs_no_cancer]
pheno_no_cancer_reg = pheno[pheno$person_id %in% IDs_no_cancer]
geno_no_cancer_reg = select.inds(geno, id %in% IDs_no_cancer)
IDs_no_cancer_reg=geno_no_cancer_reg@ped$id
length(IDs_no_cancer_reg)  # 15302
nrow(covar_no_cancer_reg) ## 15302

r_no_cancer_reg=run_PGRM_assoc(geno=geno_no_cancer_reg, pheno=pheno_no_cancer_reg,demos=covar_no_cancer_reg,covariates = c('last_age','sex','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'),
                               PGRM=PGRM, MCC=2,minimum_case_count=100,use_exclude_ranges=TRUE,LOUD=TRUE,check_sex = TRUE)
r_no_cancer_reg=annotate_results(r_no_cancer_reg,ancestry="EUR",build="hg19",calculate_power = TRUE)
get_RR(r_no_cancer_reg) ## Replicated 103 of 174 for RR=59.2%
get_AER(r_no_cancer_reg) ## Expected 694.8, replicated 479 for AE=0.689 (2766 associations for 72 uniq phecodes)

write.table(r_no_cancer_reg,file="anno_BioVU_EUR_no_cancer_reg.csv",row.names = F, col.names = T)

#####

r_cancer_only=fread(file="anno_BioVU_EUR_cancer_reg.csv",header=T,colClasses = list(character = 'phecode'))
r_no_cancer=fread(file="anno_BioVU_EUR_no_cancer_reg.csv",header=T,colClasses = list(character = 'phecode'))

r_cancer_only=r_cancer_only[,c('SNP','phecode','cases','controls','P','odds_ratio','L95','U95')]
r_no_cancer=r_no_cancer[,c('SNP','phecode','cases','controls','P','odds_ratio','L95','U95')]
r_cancer_only=annotate_results(r_cancer_only,ancestry="EUR",build="hg19",calculate_power = TRUE)
r_no_cancer=annotate_results(r_no_cancer,ancestry="EUR",build="hg19",calculate_power = TRUE)
get_RR(r_cancer_only)
get_RR(r_no_cancer)


compare_annotated_results(benchmark_results[cohort=="BioVU_EUR"],r_cancer_only)

compare_annotated_results(benchmark_results[cohort=="BioVU_EUR"],r_no_cancer)
