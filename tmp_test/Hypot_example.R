library(gaston)
library(data.table)
library(glue)
library(pgrm)


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
select.snps(geno, id %in% c('2:43698028:A:T','2:227092150:T:A'))
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

library(pgrm)
PGRM=get_PGRM(ancestry = "EUR",build="hg19")
GRS=make_GRS(PGRM,geno,phecode="244",prune=T,R2=0.2)

nrow(d)
nrow(GRS)

head(d)
head(GRS)
d=merge(d,GRS,by="person_id")

head(d)
m=glm(phecode2~GRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d)
summary(m)
summary(m)$coeff['GRS',4]

r=data.frame()

## Phecode MCC=1
#phenotype = glue('phecode (MCC=1) (n={table(d$phecode)[["1"]]})')
#phenotype
#m=glm(phecode~GRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
#summary(m)$coeff[2,4] ##  2.343958e-75
#c=confint.default(m)
#r=rbind(r,data.frame(i=1,phenotype,cases=table(d$phecode)[["1"]], controls=table(d$phecode)[["0"]],
#                     P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

## Phecode MCC=2
d$phecode2_ex = d$phecode2
d[phecode==1 & phecode2==0]$phecode2_ex=NA
phenotype = glue('phecode (MCC=2) (n={table(d$phecode2_ex)[["1"]]})')
phenotype
m=glm(phecode2_ex~GRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)
summary(m)$coeff[2,4] ##  2.343958e-75
c=confint.default(m)
r=rbind(r,data.frame(i=2,phenotype,cases=table(d$phecode2_ex)[["1"]], controls=table(d$phecode2_ex)[["0"]],
                     P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

## Meds
phenotype = glue('meds (n={table(d$meds)[["1"]]})')
phenotype
m=glm(meds~GRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)$coeff[2,4] ##  6.724116e-104
c=confint.default(m)
r=rbind(r,data.frame(i=3,phenotype,cases=table(d$meds)[["1"]], controls=table(d$meds)[["0"]],
                     P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

## labs
phenotype = glue('lab (n={table(d$lab)[["1"]]})')
phenotype
m=glm(lab~GRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)$coeff[2,4] ##  1.898473e-79
c=confint.default(m)
r=rbind(r,data.frame(i=4,phenotype,cases=table(d$lab)[["1"]], controls=table(d$lab)[["0"]],
                     P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

## PL
phenotype = glue('PL (n={table(d$PL)[["1"]]})')
phenotype
m=glm(PL~GRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)$coeff[2,4] ##  2.667879e-80
c=confint.default(m)
r=rbind(r,data.frame(i=5,phenotype,cases=table(d$PL)[["1"]], controls=table(d$PL)[["0"]],
                     P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))



## Phecode+meds+lab+PL
d$everything = 0
d[meds+lab+PL+phecode2==4]$everything = 1
d[everything==0 & PL+meds+lab+phecode>0]$everything = NA
phenotype = glue('phecode ∩ meds ∩ labs ∩ PL (n={table(d$everything)[["1"]]})')
phenotype
m=glm(everything~GRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)$coeff[2,4] ##  2.667879e-80
c=confint.default(m)
r=rbind(r,data.frame(i=6,phenotype,cases=table(d$everything)[["1"]], controls=table(d$everything)[["0"]],
                     P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

## phecode alone
#d$phecode_only = d$phecode
#d[phecode_only == 1 & meds+lab+PL>0]$phecode_only = NA
#phenotype = glue('PL (n={table(d$phecode_only)[["1"]]})')
#phenotype
#m=glm(phecode_only~GRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
#summary(m)$coeff[2,4] ##  2.667879e-80
#c=confint.default(m)
#r=rbind(r,data.frame(i=6,phenotype,cases=table(d$phecode_only)[["1"]], controls=table(d$phecode_only)[["0"]],
#                     P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

## Only phecode MCC=2
d$phecode2_only = d$phecode2_ex
d[phecode2_only == 1 & meds+lab+PL>0]$phecode2_only = NA
phenotype = glue('Only phecode (n={table(d$phecode2_only)[["1"]]})')
phenotype
m=glm(phecode2_only~GRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)$coeff[2,4] ##  2.667879e-80
c=confint.default(m)
r=rbind(r,data.frame(i=7,phenotype,cases=table(d$phecode2_only)[["1"]], controls=table(d$phecode2_only)[["0"]],
                     P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

## Only med
d$meds_only = d$meds
d[meds_only == 1 & (lab+PL+phecode2>0)]$meds_only = NA
phenotype = glue('Only med (n={table(d$meds_only)[["1"]]})')
phenotype
m=glm(meds_only~GRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)$coeff[2,4] ##  2.667879e-80
c=confint.default(m)
r=rbind(r,data.frame(i=8,phenotype,cases=table(d$meds_only)[["1"]], controls=table(d$meds_only)[["0"]],
                     P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

## Only lab
d$lab_only = d$lab
d[lab_only == 1 & meds+PL+phecode2>0]$lab_only = NA
phenotype = glue('Only lab (n={table(d$lab_only)[["1"]]})')
phenotype
m=glm(lab_only~GRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)$coeff[2,4] ##  2.667879e-80
c=confint.default(m)
r=rbind(r,data.frame(i=9,phenotype,cases=table(d$lab_only)[["1"]], controls=table(d$lab_only)[["0"]],
                     P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

## Only PL
d$PL_only = d$PL
d[PL_only == 1 & meds+lab+phecode2>0]$PL_only = NA
phenotype = glue('Only PL (n={table(d$PL_only)[["1"]]})')
phenotype
m=glm(PL_only~GRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)$coeff[2,4] ##  2.667879e-80
c=confint.default(m)
r=rbind(r,data.frame(i=10,phenotype,cases=table(d$PL_only)[["1"]], controls=table(d$PL_only)[["0"]],
                     P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

## Any
d$anything = 0
d[PL+meds+lab+phecode2>0]$anything = 1
d[anything==0 & phecode==1]$anything = NA
phenotype = glue('phecode U meds U labs ⋃ PL (n={table(d$anything)[["1"]]})')
phenotype
m=glm(anything~GRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)$coeff[2,4] ##  2.667879e-80
c=confint.default(m)
r=rbind(r,data.frame(i=10,phenotype,cases=table(d$anything)[["1"]], controls=table(d$anything)[["0"]],
                     P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))


## Phecode no med
d$phecode_no_med = d$phecode2_ex
d[!is.na(phecode_no_med) &  phecode_no_med == 1 & meds == 1]$phecode_no_med = NA
phenotype = glue('phecode ⊄ med (n={table(d$phecode_no_med)[["1"]]})')
phenotype
m=glm(phecode_no_med~GRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)$coeff[2,4] ##  2.667879e-80
c=confint.default(m)
r=rbind(r,data.frame(i=12,phenotype,cases=table(d$phecode_no_med)[["1"]], controls=table(d$phecode_no_med)[["0"]],
                     P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))


## eMERGE algorithm
d$eMERGE = 0
d[phecode2+lab+PL>0 & d$meds == 1 ]$eMERGE = 1
d[eMERGE==0 & phecode+lab+meds+PL>0]$eMERGE=NA
phenotype = glue('phecode U lab U PL ⋂ med (n={table(d$eMERGE)[["1"]]})')
phenotype
m=glm(eMERGE~GRS+sex+last_age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8,data=d,family="binomial")
summary(m)$coeff[2,4] ##  2.667879e-80
c=confint.default(m)
r=rbind(r,data.frame(i=13,phenotype,cases=table(d$eMERGE)[["1"]], controls=table(d$eMERGE)[["0"]],
                     P=summary(m)$coeff[2,4] , Beta=summary(m)$coeff[2,1], L95=c[2,1], U95=c[2,2]))

r


r$P_lab=glue('P={formatC(r$P, format = "e",digits=1)}')

ggplot(r, aes(x=exp(Beta),y=reorder(phenotype,i),label=P_lab))+geom_errorbar(aes(xmin=exp(L95), xmax=exp(U95)), width=.1)+
  geom_point()+geom_text(hjust=-.4, vjust=-1)+ylab('')+xlab('GRS Beta')+scale_y_discrete(limits=rev)+theme_classic()+
  geom_vline(xintercept=1,color="grey",linetype="dashed")


