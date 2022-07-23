library(gaston)
library(data.table)

geno_file="~/Dropbox (VUMC)/PGRM/code/data/PGRM_MEGA4_Eur_adult_V1"
covar_file = "~/Dropbox (VUMC)/PGRM/code/data/covariates_MEGA4_EUR.csv"
pheno_file = "~/Dropbox (VUMC)/PGRM/code/data/phecodes_V1_2_UP.csv"
hypot_FT4 = "~/Dropbox (VUMC)/PGRM/code/data/MEGA_low_FT4.csv"
hypot_TSH = "~/Dropbox (VUMC)/PGRM/code/data/MEGA_high_TSH.csv"
hypot_PL = "~/Dropbox (VUMC)/PGRM/code/data/MEGA_hypot_PL.csv"
hypot_meds = "~/Dropbox (VUMC)/PGRM/code/data/MEGA_hypot_med_2x.csv"

## Read genotype file
geno=read.bed.matrix(geno_file)
geno=select.snps(geno, maf > 0.01)
geno=select.snps(geno, callrate > 0.98)
IDs=geno@ped$id

## Read covariate file
covar=read.csv(covar_file,header=TRUE,stringsAsFactors = F)
names(covar)[1]="person_id"
covar=data.table(covar,key="person_id")

ages=read.csv("~/Dropbox (VUMC)/PGRM/code/data/covariates_MEGA4_ages.csv",header=TRUE,stringsAsFactors = F)
names(ages)[1]="person_id"
ages=data.table(ages,key="person_id")
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

## read labs

TSH = fread(hypot_TSH)
setkeyv(TSH, c( 'person_id'))
TSH$TSH=1

FT4 = fread(hypot_FT4)
setkeyv(FT4, c( 'person_id'))
FT4$FT4=1

PL = fread(hypot_PL)
setkeyv(PL, c( 'person_id'))
PL$PL = 1

meds = fread(hypot_meds)
setkeyv(meds, c( 'person_id'))
meds$meds = 1

phecode = pheno[phecode=="244",c("person_id","N")]
names(phecode)[2]="phecode_count"
phecode$phecode = 1

phecode2 = pheno[phecode=="244" & N>=2,c("person_id")]
phecode2$phecode2 = 1



d = merge(covar, TSH, by="person_id", all.x=T)
d = merge(d, FT4, by="person_id", all.x=T)
d = merge(d, PL, by="person_id", all.x=T)
d = merge(d, meds, by="person_id", all.x=T)
d = merge(d, phecode, by="person_id", all.x=T)
d = merge(d, phecode2, by="person_id", all.x=T)
head(d)

d[is.na(TSH)]$TSH= 0
d[is.na(FT4)]$FT4= 0
d[is.na(PL)]$PL= 0
d[is.na(meds)]$meds= 0
d[is.na(phecode)]$phecode= 0
d[is.na(phecode2)]$phecode2= 0
head(d)

d$lab = 0
d[TSH == 1 | FT4==1]$lab = 1

table(d$lab)
table(d$meds)
table(d$phecode)
table(d$phecode2)
table(d$PL)

x = list(
  phecode = d[phecode2==1,]$person_id,
  meds = d[meds==1,]$person_id,
  lab = d[lab==1,]$person_id,
  PL = d[PL==1,]$person_id
)

library(ggvenn)
ggvenn(x,  stroke_size = 0.5, show_percentage = T )

source('~/Dropbox (VUMC)/PGRM/code/run_PGRM_functions.R')
PGRM=get_PGRM(ancestry = "EUR",build="hg19")

PRS=make_PRS(PGRM,geno,"244",prune=T,R2=0.2)






#d$any_meas = 0
#d[TSH_meas == 1 | FT4_meas==1]$any_meas = 1
#table(d$any_meas)


## read meds

#med=read.csv("data/MEGA_hypot_med_lim.csv",header=T)
med=read.csv("data/MEGA_hypot_med_2x.csv",header=T)
nrow(med)
med$med = 1
d = merge(d, med, by="ID", all.x=T)
d[is.na(med)]$med= 0
table(d$med)
table(d$med,d$any_lab)
table(d$med,d$any_meas)

## set pheno
phecode=pheno[phecode=="244" & N>=1,'ID']
phecode$phecode = 1
d = merge(d, phecode, by="ID", all.x=T)
d[is.na(phecode)]$phecode= 0

## set pheno_MCC2
phecode2=pheno[phecode=="244" & N>=2,'ID']
phecode2$phecode_MCC2 = 1
d = merge(d, phecode2, by="ID", all.x=T)
d[is.na(phecode_MCC2)]$phecode_MCC2= 0
d[phecode_MCC2 == 0 & phecode==1]$phecode_MCC2=NA
table(d$phecode,d$phecode_MCC2)

## set PRS
PRS=make_PRS(PGRM,geno,"244",prune=T,R2=0.2)
d=merge(d,PRS[,c('ID','PRS')],by="ID")

x = list(
  phecode = d[phecode==1,]$ID,
  med = d[med==1,]$ID,
  TSH = d[TSH==1,]$ID,
  FT4 = d[FT4==1,]$ID
)

#x = list(
#  phecode = d[phecode==1,]$ID,
#  med = d[med==1,]$ID,
#  any_lab = d[any_lab==1,]$ID
#)


library(ggvenn)
ggvenn(x,  stroke_size = 0.5, show_percentage = T )

d$phe_OR_lab_plus_med = 0
d[(any_lab==1 | phecode == 1) & med==1]$phe_OR_lab_plus_med = 1
table(d$phe_OR_lab_plus_med,useNA = "always")

d$phe2_OR_lab_plus_med = 0
d[(any_lab==1 | phecode_MCC2 == 1) & med==1]$phe2_OR_lab_plus_med = 1
d[phe2_OR_lab_plus_med == 0 & (any_lab==1 | phecode == 1 | med == 1)]$phe2_OR_lab_plus_med =NA
table(d$phe2_OR_lab_plus_med,useNA = "always")

d$gold_all = 0
d[any_lab==1 & phecode == 1 & med==1]$gold_all = 1
d[gold_all == 0 & (any_lab==1 | phecode == 1 | med == 1)]$gold_all =NA
table(d$gold_all,useNA = "always")

## test phecode MCC2

r=data.frame()

table(d$phecode)
phenotype = "Phecode, MCC=1 (n=11,961)"
m=glm(phecode~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)
summary(m)$coeff[2,4] ##  3.206464e-72
c=confint.default(m)
r=rbind(r,data.frame(i=1,phenotype,cases=table(d$phecode)[["1"]], controls=table(d$phecode)[["0"]], P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

table(d$phecode_MCC2)
phenotype = "Phecode, MCC=2 (n=9,068)"
m=glm(phecode_MCC2~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)
summary(m)$coeff[2,4] ##  2.708795e-78
c=confint.default(m)
r=rbind(r,data.frame(i=2,phenotype,cases=table(d$phecode)[["1"]], controls=table(d$phecode)[["0"]], P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

table(d$med)
phenotype = "Medication use (n=12,782)"
m=glm(med~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)
summary(m)$coeff[2,4] ##  2.708795e-78
c=confint.default(m)
r=rbind(r,data.frame(i=3,phenotype,cases=table(d$phecode)[["1"]], controls=table(d$phecode)[["0"]], P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))


table(d$TSH)
phenotype = "TSH high (n=6,873)"
m=glm(TSH~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)
summary(m)$coeff[2,4] ##  2.708795e-78
c=confint.default(m)
r=rbind(r,data.frame(i=4,phenotype,cases=table(d$phecode)[["1"]], controls=table(d$phecode)[["0"]], P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

table(d$FT4)
phenotype = "Low FT4 (n=784)"
m=glm(FT4~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)
summary(m)$coeff[2,4] ##  2.708795e-78
c=confint.default(m)
r=rbind(r,data.frame(i=5,phenotype,cases=table(d$phecode)[["1"]], controls=table(d$phecode)[["0"]], P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

table(d$any_lab)
phenotype = "Abnormal lab (n=7089)"
m=glm(any_lab~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)
summary(m)$coeff[2,4] ## 3.597616e-81
c=confint.default(m)
r=rbind(r,data.frame(i=6,phenotype,cases=table(d$phecode)[["1"]], controls=table(d$phecode)[["0"]], P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))



## test gold
table(d$phe_OR_lab_plus_med,useNA="always")  # 9231
phenotype = "(Phecode | lab) + med (n=10,332)"
m=glm(phe_OR_lab_plus_med~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)
summary(m)$coeff[2,4] ##  3.733946e-98
c=confint.default(m)
r=rbind(r,data.frame(i=7,phenotype,cases=table(d$phecode)[["1"]], controls=table(d$phecode)[["0"]], P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

table(d$phe2_OR_lab_plus_med,useNA="always")  # 9231
phenotype = "Phecode2 | lab + (med) (n=9096)"
m=glm(phe_OR_lab_plus_med~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)
summary(m)$coeff[2,4] ##  3.733946e-98
c=confint.default(m)
r=rbind(r,data.frame(i=8,phenotype,cases=table(d$phecode)[["1"]], controls=table(d$phecode)[["0"]], P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

table(d$gold_all,useNA="always")
phenotype = "Phecode + med + lab (n=4666)"
m=glm(gold_all~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)
summary(m)$coeff[2,4] ##  3.733946e-98
c=confint.default(m)
r=rbind(r,data.frame(i=9,phenotype,cases=table(d$phecode)[["1"]], controls=table(d$phecode)[["0"]], P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

ggplot(r, aes(x=exp(Beta),y=reorder(phenotype,i)))+geom_errorbar(aes(xmin=exp(L95), xmax=exp(U95)), width=.1)+geom_point()+ylab('')+scale_y_discrete(limits=rev)+theme_classic()+geom_vline(xintercept=1,color="grey",linetype="dashed")


rc[2,]
#    2.5 %    97.5 %
# 0.8810339 1.0621652

## Rand 5%
cases = d[gold==1]
cases_ID= cases[sample(nrow(cases)*0.05)]$ID
controls=d[gold==0]
control_ID=controls[sample(length(cases_ID))]$ID
d$gold_rand05 = d$gold
d[ID %in% cases_ID]$gold_rand05=0
d[ID %in% control_ID]$gold_rand05=1
table(d$gold,d$gold_rand05)
m=glm(gold_rand05~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)
summary(m)$coeff[2,4] ##  1.911651e-85
c=confint.default(m)
c[2,]
#     2.5 %     97.5 %
# 0.8101696 0.9903139

run_rand_assoc= function(d,rand=0,swap=F){
  len=round(nrow(d[gold==1])*rand,0)

  total_number_of_cases = nrow(d[gold==1])
  cases_ID= sample(d[gold==1]$ID,size=len,replace=F)
  control_ID_total=sample(d[gold==0]$ID,size=total_number_of_cases,replace=F)
  control_ID=sample(control_ID_total,size=len,replace=F)


  d$gold_rand = d$gold
  if(swap==T){
    d[ID %in% cases_ID]$gold_rand=0
  } else {
    d[ID %in% cases_ID]$gold_rand=NA
    d[gold == 0 & ID %in% control_ID_total]$gold_rand = NA
  }
  d[ID %in% control_ID]$gold_rand=1
  table(d$gold,d$gold_rand,useNA = "always")
  m=glm(gold_rand~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
  summary(m)
  P=summary(m)$coeff[2,4] ##  2.160914e-77
  Beta=summary(m)$coeff[2,1]
  c=confint.default(m)
  L95=c[2,1]
  U95=c[2,2]
  Cstat=Cstat(m)
  results = data.frame(Rand=rand,Perecent_random= paste0(rand*100, "%"),
                       cases=nrow(d[gold_rand==1]), controls=nrow(d[gold_rand==0]), P=P, Beta = Beta, L95=L95, U95=U95,Cstat=Cstat)
  return(results)
}

randomize_pheno = function(d,rand=0,swap=T){
  len=round(nrow(d[gold==1])*rand,0)

  total_number_of_cases = nrow(d[gold==1])
  cases_ID= sample(d[gold==1]$ID,size=len,replace=F)
  control_ID_total=sample(d[gold==0]$ID,size=total_number_of_cases,replace=F)
  control_ID=sample(control_ID_total,size=len,replace=F)
  d$gold_rand = d$gold
  if(swap==T){
    d[ID %in% cases_ID]$gold_rand=0
  } else {
    d[ID %in% cases_ID]$gold_rand=NA
    d[gold == 0 & ID %in% control_ID_total]$gold_rand = NA
  }
  d[ID %in% control_ID]$gold_rand=1
  return(d)
}

#func = function(DataSet,indices,rand=0.05){
#  d=DataSet[indices,]
#  d=randomize_pheno(d,rand=0.1,swap=T)
#  fit = glm(gold_rand~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
#  return(DescTools::Cstat(fit))
#}
#library(boot)
#library(DescTools)
#res = boot(d,statistic=func,R=1000)
#plot(res)
#boot.ci(res,type="norm")


foo=randomize_pheno(d,rand=0.1,swap=T)

run_rand_assoc(d,rand=0.1,swap=F)
run_rand_assoc(d,rand=0.1,swap=T)

rands=seq(from=0, to=1,by=.05)
r_rand = data.frame()
for(i in 1:length(rands)){
  r_rand=rbind(r_rand,run_rand_assoc(d,rand=rands[i],swap=F))
}
r_rand
write.table(r_rand,"results/BioVU_EUR_hypot_rand.txt",sep="\t",row.names = F,col.names = T)








library(ggplot2)
r_rand$Rand=as.factor(r_rand$Rand)
ggplot(r_rand, aes(x=Beta, y=Rand)) +geom_errorbar(aes(xmin=L95, xmax=U95), width=.1)+geom_point()+scale_y_discrete(limits=rev)+theme_bw()




d1=randomize_pheno(d,rand=0,swap=F)
d2=randomize_pheno(d,rand=0.05,swap=F)
m1=glm(gold_rand~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d1,family="binomial")
m2=glm(gold_rand~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d2,family="binomial")

d1$set = "rand_00"
d2$set = "rand_01"
d_full=rbind(d1,d2)
m=glm(gold_rand~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+set,data=d_full,family="binomial")
summary(m)

AIC(m1,m2)
anova(m1, m2, test = "LRT")
Cstat(m1)
Cstat(m2)


func = function(DataSet,indices){
  d=DataSet[indices,]
  fit = glm(gold_rand~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
  return(DescTools::Cstat(fit))
}
res = boot(d1,statistic=func,R=100,parallel="multicore")
plot(res)
boot.ci(res,type="norm")
#95%   ( 0.6890,  0.7003 )

d2=randomize_pheno(d,rand=0.05,swap=F)
res2 = boot(d2,statistic=func,R=100)
plot(res2)
boot.ci(res2,type="norm")
# 95%   ( 0.6783,  0.6884 )

GetBootStrap_Pvalue <- function(formula,DataSet1,DataSet2,indices) {
  d1 <- DataSet1[indices,]
  d2 <- DataSet2[indices,]
  fit1 = glm(formula,data=d1, family="binomial")
  result1=DescTools::Cstat(fit1)
  fit2 = glm(formula,data=d2, family="binomial")
  result2=DescTools::Cstat(fit2)
  return(result1-result2)
}
formula=as.formula("gold_rand~PRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8")
