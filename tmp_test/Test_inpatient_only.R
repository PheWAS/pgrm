library(data.table)
library(gaston)
library(pgrm)
library(glue)
geno_file="~/Dropbox (VUMC)/PGRM/code/data/PGRM_MEGA4_Eur_adult_V1"
covar_file = "~/Dropbox (VUMC)/PGRM/code/data/covariates_MEGA4_EUR.csv"
pheno_file = "~/Dropbox (VUMC)/PGRM/code/data/phecodes_V1_2_UP_I.csv"

library(pgrm)
results_BBJ
results_BioVU_AFR
results_BioVU_EUR
results_MGI
results_UKBB


## Load PGRM
PGRM=get_PGRM(ancestry="EUR",build="hg19")

cohort_info = read.csv(file="~/Dropbox (VUMC)/PGRM/code/PGRM/PGRM_cohort.csv",stringsAsFactors = F)

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

r=run_PGRM_assoc(geno=geno, pheno=pheno,demos=covar,covariates = c('last_age','sex','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'),
                 PGRM=PGRM, MCC=2,minimum_case_count=100,use_exclude_ranges=TRUE,LOUD=TRUE,check_sex = TRUE)

r_anno=annotate_results(r,ancestry="EUR",build="hg19",calculate_power = TRUE)
r_anno=merge(r_anno,cohort_info,by="assoc_ID")
get_RR(r_anno[BioVU_only==0])
## Replicated 280 of 365 for RR=76.7%
get_AER(r_anno[BioVU_only==0])
## Expected 995.1, replicated 828 for AE=0.832 (2804 associations for 71 uniq phecodes)
nrow(r_anno[Power>0.8])
nrow(r_anno[cases>=cases_needed])

new_BioVU_EUR_INPT=r_anno[,c('SNP','phecode','cases','controls','odds_ratio','P','L95','U95','BioVU_only')]
names(new_BioVU_EUR_INPT)[9]="cohort_match"
write.table(new_BioVU_EUR_INPT,file="results_BioVU_EUR_INPT.csv",row.names = F, col.names = T)

