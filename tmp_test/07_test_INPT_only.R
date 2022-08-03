library(data.table)
library(gaston)
library(pgrm)
library(glue)
geno_file="~/Dropbox (VUMC)/PGRM/code/data/PGRM_MEGA4_Eur_adult_V1"
covar_file = "~/Dropbox (VUMC)/PGRM/code/data/covariates_MEGA4_EUR.csv"
pheno_file = "~/Dropbox (VUMC)/PGRM/code/data/phecodes_V1_2_UP_I.csv"

set.seed(4)

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


foo=annotate_results(results_BioVU_EUR,ancestry="EUR")
PGRM=PGRM[!assoc_ID %in% foo[cohort_match==1]$assoc_ID]
get_RR(foo)

r_inpt_only=run_PGRM_assoc(geno=geno, pheno=pheno,demos=covar,covariates = c('last_age','sex','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'),
                 PGRM=PGRM, MCC=2,minimum_case_count=100,use_exclude_ranges=TRUE,LOUD=TRUE,check_sex = TRUE)
r_inpt_only=annotate_results(r_inpt_only,ancestry="EUR",build="hg19",calculate_power = TRUE)
get_RR(r_inpt_only) ## Replicated 286 of 375 for RR=76.3%
get_AER(r_inpt_only) ## Expected 1014, replicated 828 for AE=0.817 (2805 associations for 71 uniq phecodes)

write.table(r_inpt_only,file="anno_BioVU_EUR_INPT_only.csv",row.names = F, col.names = T)


#######

r_inpt_only=fread(file="anno_BioVU_EUR_INPT_only.csv",header=T,colClasses = list(character = 'phecode'))

r_inpt_only=r_inpt_only[,c('SNP','phecode','cases','controls','P','odds_ratio','L95','U95')]

r_inpt_only=annotate_results(r_inpt_only,ancestry="EUR",build="hg19",calculate_power = TRUE)
nrow(r_inpt_only)
get_RR(r_inpt_only,include="all") ##
get_powered_rate(r_inpt_only)
get_RR(r_inpt_only,include="powered") ##
get_AER(r_inpt_only) ##
r_inpt_only[is.na(Power)]

prop.table(table(r_inpt_only$CI_overlap))

t.test(r_inpt_only$cat_OR,r_inpt_only$rOR,paired=T) ## mean diff 0.08509277
t.test(r_inpt_only[rep==1]$cat_OR,r_inpt_only[rep==1]$rOR,paired=T) ## mean diff 0.06459039


compare_annotated_results(benchmark_results[cohort=="BioVU_EUR"],r_inpt_only)
