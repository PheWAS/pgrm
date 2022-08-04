library(data.table)
library(gaston)
library(pgrm)
library(glue)
geno_file="~/Dropbox (VUMC)/PGRM/code/data/PGRM_MEGA4_Eur_adult_V1"
covar_file = "~/Dropbox (VUMC)/PGRM/code/data/covariates_MEGA4_EUR.csv"
pheno_file = "~/Dropbox (VUMC)/PGRM/code/data/phecodes_V1_2_UP.csv"
med_home_file = "~/Dropbox (VUMC)/PGRM/code/data/MEGA_med_home.csv"

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


## read med home info
med_home = read.csv(med_home_file,header=TRUE,stringsAsFactors = F)
med_home = data.table(med_home)
med_home = med_home[person_id %in% IDs]
med_home = med_home[person_id %in% covar$person_id]
nrow(med_home) ## 59373
nrow(med_home[cpt_date>=2]) # 50189
nrow(med_home[cpt_year>=2]) # 44712
nrow(med_home[cpt_date>=3]) # 32594
nrow(med_home[cpt_date>=4]) # 24497

med_home_IDs = sample(med_home[cpt_year>=4]$person_id)
length(med_home_IDs)
nrow(covar[!person_id %in% med_home_IDs])


## read phenotype file
pheno=read.csv(pheno_file,header=TRUE,colClasses=c("character","character","integer"),stringsAsFactors = F)
names(pheno)[1]="person_id"
names(pheno)[2]="phecode"
pheno=pheno[pheno$person_id %in% med_home_IDs,] ## include only individuals with genotype data
pheno=data.table(pheno,key="person_id")
nrow(pheno)



covar_med_home = covar[covar$person_id %in% med_home_IDs]
pheno_med_home = pheno[pheno$person_id %in% med_home_IDs]
geno_med_home = select.inds(geno, id %in% med_home_IDs)
IDs_med_home=geno_med_home@ped$id
length(IDs_med_home)  # 44712
nrow(covar_med_home) ## 44712



r_med_home=run_PGRM_assoc(geno=geno_med_home, pheno=pheno_med_home,demos=covar_med_home,covariates = c('last_age','sex','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'),
                 PGRM=PGRM, MCC=2,minimum_case_count=100,use_exclude_ranges=TRUE,LOUD=TRUE,check_sex = TRUE)
r_med_home=annotate_results(r_med_home,ancestry="EUR",build="hg19",calculate_power = TRUE)
get_RR(r_med_home) ## Replicated 514 of 710 for RR=72.4%
get_AER(r_med_home) ## Expected 1525, replicated 1205 for AE=0.79 (3255 associations for 103 uniq phecodes)

write.table(r_med_home,file="anno_BioVU_EUR_med_home4.csv",row.names = F, col.names = T)


#####

r_med_home4=fread(file="anno_BioVU_EUR_med_home4.csv",header=T,colClasses = list(character = 'phecode'))


r_med_home4=r_med_home4[,c('SNP','phecode','cases','controls','P','odds_ratio','L95','U95')]
r_med_home4=annotate_results(r_med_home4,ancestry="EUR",build="hg19",calculate_power = TRUE)



r_med_home2=fread(file="anno_BioVU_EUR_med_home2.csv",header=T,colClasses = list(character = 'phecode'))


r_med_home2=r_med_home2[,c('SNP','phecode','cases','controls','P','odds_ratio','L95','U95')]
r_med_home2=annotate_results(r_med_home2,ancestry="EUR",build="hg19",calculate_power = TRUE)

get_RR(r_new)
get_RR(r_new2)

compare_annotated_results(benchmark_results[cohort=="BioVU_EUR"],r_med_home2)

compare_annotated_results(benchmark_results[cohort=="BioVU_EUR"],r_med_home4)

t=t.test(r$rOR.x,r$rOR.y,paired=T)


r1=benchmark_results[cohort=="BioVU_EUR",c('assoc_ID','rsID','phecode','phecode_string','category_string','odds_ratio','P','L95','U95','rep','powered','Power','rOR','rL95','rU95')]
r2=r_med_home[,c('assoc_ID','odds_ratio','P','L95','U95','rep','powered','Power','rOR','rL95','rU95')]
r=merge(r1,r2,by="assoc_ID")
head(r)
mean(r$rOR.x)
mean(r$rOR.y)

r$both_rep = 0
r[rep.x==1 & rep.y==1]$both_rep = 1

mean(r[both_rep==1]$rOR.x)-mean(r[both_rep==1]$rOR.y)
mean(r[both_rep==1]$rOR.y)
t.test(r[both_rep==1]$rOR.x,r[both_rep==1]$rOR.y,paired=TRUE)
