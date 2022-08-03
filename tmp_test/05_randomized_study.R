library(data.table)
library(gaston)
library(pgrm)
library(glue)
geno_file="~/Dropbox (VUMC)/PGRM/code/data/PGRM_MEGA4_Eur_adult_V1"
covar_file = "~/Dropbox (VUMC)/PGRM/code/data/covariates_MEGA4_EUR.csv"
pheno_file = "~/Dropbox (VUMC)/PGRM/code/data/phecodes_V1_2_UP.csv"

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

#sex_specific=c('182','185','615','618','180.1','256.4','187.2','626.12','174.11','174.2','184.11','218.1')
#PGRM=PGRM[!phecode %in% sex_specific]

foo=annotate_results(results_BioVU_EUR,ancestry="EUR")
PGRM=PGRM[!assoc_ID %in% foo[cohort_match==1]$assoc_ID]

#PGRM=head(PGRM,n=20)


run_rand_assoc = function(rand_frac){
  ID_map = data.table(sample(covar$person_id))
  names(ID_map)[1]="person_id"
  ID_len=nrow(ID_map)
  ID_map$person_id_rand = ID_map$person_id
  ID_map[1:round(ID_len*rand_frac,0),]$person_id_rand=sample(ID_map[1:round(ID_len*rand_frac,0),]$person_id_rand,replace=F)
  nrow(ID_map[person_id != person_id_rand])
  nrow(ID_map)
  setkeyv(ID_map, c('person_id'))

  c = merge(covar, ID_map,by="person_id")
  c$person_id=c$person_id_rand

  p = merge(pheno, ID_map,by="person_id")
  p$person_id=p$person_id_rand
  p=p[,1:3]

  r=run_PGRM_assoc(geno=geno, pheno=p,demos=c,covariates = c('last_age','sex','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'),
                 PGRM=PGRM, MCC=2,minimum_case_count=100,use_exclude_ranges=TRUE,LOUD=TRUE,check_sex = TRUE)
  r=annotate_results(r,ancestry="EUR",build="hg19",calculate_power = TRUE)
  return(r)
}


#r000=run_rand_assoc(rand_frac=0)
#get_RR(r000[!assoc_ID %in% foo[cohort_match==1]$assoc_ID])

r100=run_rand_assoc(rand_frac=1.0)
write.table(r100,file="tmp_test/OUT/r100.csv",row.names = F, col.names = T)


r075=run_rand_assoc(rand_frac=0.75)
write.table(r075,file="tmp_test/OUT/r075.csv",row.names = F, col.names = T)



r050=run_rand_assoc(rand_frac=0.50)
write.table(r050,file="tmp_test/OUT/r050.csv",row.names = F, col.names = T)


r025=run_rand_assoc(rand_frac=0.25)
write.table(r025,file="tmp_test/OUT/r025.csv",row.names = F, col.names = T)


r000=run_rand_assoc(rand_frac=0)
write.table(r000,file="tmp_test/OUT/r000.csv",row.names = F, col.names = T)


r010=run_rand_assoc(rand_frac=0.10)
write.table(r010,file="tmp_test/OUT/r010.csv",row.names = F, col.names = T)


r020=run_rand_assoc(rand_frac=0.20)
write.table(r020,file="tmp_test/OUT/r020.csv",row.names = F, col.names = T)


r030=run_rand_assoc(rand_frac=0.30)
write.table(r030,file="tmp_test/OUT/r030.csv",row.names = F, col.names = T)


r040=run_rand_assoc(rand_frac=0.40)
write.table(r040,file="tmp_test/OUT/r040.csv",row.names = F, col.names = T)



get_RR(r100)
get_AER(r100)

get_RR(r075)
get_AER(r075)

get_RR(r050)
get_AER(r050)

get_RR(r025)
get_AER(r025)

get_RR(r000)
get_AER(r000)

get_RR(r010)
get_AER(r010)

get_RR(r020)
get_AER(r020)

get_RR(r030)
get_AER(r030)

get_RR(r040)
get_AER(r040)
