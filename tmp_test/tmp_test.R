library(data.table)
library(gaston)
library(pgrm)
library(glue)
#geno_file="~/Dropbox (VUMC)/PGRM/code/data/PGRM_MEGA4_Eur_adult_V1"
#covar_file = "~/Dropbox (VUMC)/PGRM/code/data/covariates_MEGA4_EUR.csv"
geno_file="~/Dropbox (VUMC)/PGRM/code/data/PGRM_MEGA4_Afr_adult_V1"
covar_file = "~/Dropbox (VUMC)/PGRM/code/data/covariates_MEGA4_AFR.csv"
pheno_file = "~/Dropbox (VUMC)/PGRM/code/data/phecodes_V1_2_UP.csv"


## Load PGRM
PGRM=get_PGRM(ancestry="AFR",build="hg19")

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

## read phenotype file
pheno=read.csv(pheno_file,header=TRUE,colClasses=c("character","character","integer"),stringsAsFactors = F)
names(pheno)[1]="person_id"
names(pheno)[2]="phecode"
pheno=pheno[pheno$person_id %in% IDs,] ## include only individuals with genotype data
pheno=data.table(pheno,key="person_id")

demos=covar

MCC=2
minimum_case_count=100
use_exclude_ranges=TRUE
LOUD=TRUE
check_sex=TRUE

r=run_PGRM_assoc(geno=geno, pheno=pheno,demos=demos,covariates = c('last_age','sex','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'),
               PGRM=PGRM, MCC=2,minimum_case_count=100,use_exclude_ranges=TRUE,LOUD=TRUE,check_sex = TRUE)

r_anno=annotate_results(r,ancestry="AFR",build="hg19",calculate_power = TRUE)
nrow(r_anno)
head(r_anno)
r_anno=merge(r_anno,cohort_info,by="assoc_ID")
get_RR(r_anno[BioVU_only==0])
get_AER(r_anno[BioVU_only==0])
r_anno[Power>0.8]
table(r_anno[Power>0.8]$rep)
r_anno[powered>0.8 & rep == 0]

head(results_BioVU_AFR)
head(r_anno)
cur_SNP="7:15434230:G:C"
r_anno[SNP==cur_SNP]
results_BioVU_AFR[SNP==cur_SNP]
r_anno2=annotate_results(results_BioVU_AFR,ancestry="AFR",build="hg19",calculate_power = TRUE)
get_RR(r_anno2)

UKBB_anno=annotate_results(results_UKBB[cohort_match==0],ancestry="EUR",build="hg38",calculate_power = TRUE)
MGI_anno=annotate_results(results_MGI[cohort_match==0],ancestry="EUR",build="hg38",calculate_power = TRUE)
BBJ_anno=annotate_results(results_BBJ[cohort_match==0],ancestry="EAS",build="hg19",calculate_power = TRUE)
get_RR(UKBB_anno)
get_RR(MGI_anno)
get_RR(BBJ_anno)
prop.table(table(UKBB_anno[Power>0.8]$rep))
prop.table(table(MGI_anno[Power>0.8]$rep))
prop.table(table(BBJ_anno[Power>0.8]$rep))
table(BBJ_anno[Power>0.8]$rep)

new_BioVU_AFR=r_anno[,c('SNP','phecode','cases','controls','odds_ratio','P','L95','U95','BioVU_only')]
names(new_BioVU_AFR)[9]="cohort_match"
write.table(results_BioVU_EUR,file="results_BioVU_EUR.csv",row.names = F, col.names = T)


r_anno2=annotate_results(results_BioVU_EUR[cohort_match==0],ancestry="EUR",build="hg19",calculate_power = TRUE)
head(r_anno2)
get_RR(r_anno2)

## 387
r
PGRM_ALL[phecode=="401"]
head(r)


#run_PGRM_assoc = function(geno, pheno, demos,covariates, PGRM,MCC=2,minimum_case_count=100,use_exclude_ranges=TRUE,check_sex=FALSE,LOUD=TRUE){

sex_check_phecode = function(phecode){
  if(phecode %in% pgrm::phecode_info[sex=="Male"]$phecode){
    return("M")
  }
  if (phecode %in% pgrm::phecode_info[sex=="Female"]$phecode){
    return("F")
  }
  return("B")}

covariates = c('last_age','sex','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8')

  # ## create formulas for glm using covariates; formula_string_no_sex is for sex-specific phenotypes
covar_list = paste(covariates, collapse ='+')
formula_string=paste("pheno~genotype+", covar_list, sep='')
formula_string_no_sex=gsub("\\+sex","",formula_string)
formula_string_no_sex=gsub("sex\\+","",formula_string_no_sex)
if(LOUD==TRUE){
  print(glue('Formula: {formula_string}'))}
formula_string=as.formula(formula_string)
formula_string_no_sex=as.formula(formula_string_no_sex)

## Change person_id to character to ensure merging works
demos$person_id=as.character(demos$person_id)
pheno$person_id=as.character(pheno$person_id)

# filter demos, pheno, and geno for intersection of IDs
IDs=union(geno@ped$id,demos$person_id)
demos =demos[person_id %in% IDs]
#pheno=pheno[person_id %in% IDs]
geno=select.inds(geno, id %in% IDs)
pheno =pheno[person_id %in% IDs]

## get available phecodes list
phecode_counts=pheno[N>=MCC, .(cases = .N), by=phecode]
available_phecodes = phecode_counts[cases>=minimum_case_count]$phecode

# filter PGRM for available SNPs and phecodes
SNPs=geno@snps$id
PGRM=PGRM[PGRM$SNP %in% SNPs & PGRM$phecode %in% available_phecodes,]

assoc_to_run=unique(PGRM[,c("SNP", "phecode")])
assoc_num = nrow(assoc_to_run)
if(assoc_num == 0){
  print("There are no eligable associations to run")
  return(0)}
if(LOUD==TRUE){
  print(glue('Attempting {assoc_num} association tests'))
}

results=data.frame()

  for(i in 1:nrow(assoc_to_run)){
    cur_SNP = assoc_to_run[i,]$SNP
    cur_phecode = assoc_to_run[i,]$phecode
    cur_SNP_index=which(SNPs %in% c(cur_SNP))

    if(LOUD==TRUE){
      print(glue('{i} SNP: {cur_SNP} Phecode: {cur_phecode}'))
    }


    g=data.table(as.matrix(geno[,cur_SNP_index]),keep.rownames=TRUE)
    names(g)[1]='person_id'
    names(g)[2]="genotype"
    setkey(g,"person_id")

    g$genotype = abs(g$genotype-2) ## gaston codes things "backwards" from plink. 2== HOM for ref. Flip this around
    d$person_id=as.character(demos$person_id)
    d=merge(g,demos,by="person_id")

    d=get_pheno(pheno=pheno,demos=d, phecode=cur_phecode,MCC=MCC,
                use_exclude_ranges=use_exclude_ranges, check_sex=check_sex)

    if(check_sex==TRUE){
      phecode_sex=sex_check_phecode(cur_phecode)
      print(phecode_sex)
      if(phecode_sex !="B"){
        formula_string=formula_string_no_sex
      }
    }

    n_case=nrow(d[pheno==1,])
    n_control=nrow(d[pheno==0,])
    if(n_case < minimum_case_count ) {
      next}

    m=glm(formula_string, data=d,family="binomial")
    #  conf=confint.default(m)
    #  P=summary(m)$coeff[2,4]
    #  odds_ratio=exp(summary(m)$coeff[2,1])

    #  L95=exp(conf[2,1])
    #  U95=exp(conf[2,2])
    #  result = data.frame(SNP=cur_SNP, phecode=cur_phecode,cases=n_case,controls=n_control, P=P, odds_ratio=odds_ratio, L95=L95, U95=U95)
    #  results = rbind(results, result)
  }
  return(results)
}

