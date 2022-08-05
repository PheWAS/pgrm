##foobar
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
pheno=data.table(pheno,key="person_id")
nrow(pheno)


for(i in 2:2){
  outfile=glue("~/Documents/GitHub/pgrm/tmp_test/OUT/results_BioVU_EUR_{i}.csv")
  r=run_PGRM_assoc(geno=geno, pheno=pheno,demos=covar,covariates = c('last_age','sex','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8'),PGRM=PGRM, MCC=i,minimum_case_count=100,use_exclude_ranges=TRUE,LOUD=TRUE,check_sex = TRUE)
  r=annotate_results(r,ancestry="EUR",build="hg19",calculate_power = TRUE)
  write.table(r,file=outfile,row.names = F, col.names = T)
}



r_MCC1=r
get_RR(r_MCC1,include="all")
get_RR(r_MCC1)

r_MCC=data.table()
for(i in 1:8){
  infile=glue("~/Documents/GitHub/pgrm/tmp_test/OUT/results_BioVU_EUR_{i}.csv")
  r=fread(infile)
  r$MCC= i
  r_MCC=rbind(r_MCC,r)
  #write.table(r,file=outfile,row.names = F, col.names = T)
}
MCC_summary=data.table(MCC=1, RR=get_RR(r_MCC[MCC==1]), RR_all=get_RR(r_MCC[MCC==1], include="all"), Power=get_powered_rate(r_MCC[MCC==1]),AER=get_AER(r_MCC[MCC==1]))
MCC_summary=rbind(MCC_summary,data.table(MCC=2, RR=get_RR(r_MCC[MCC==2]), RR_all=get_RR(r_MCC[MCC==2], include="all"), Power=get_powered_rate(r_MCC[MCC==2]), AER=get_AER(r_MCC[MCC==2])))
MCC_summary=rbind(MCC_summary,data.table(MCC=3, RR=get_RR(r_MCC[MCC==3]), RR_all=get_RR(r_MCC[MCC==3], include="all"), Power=get_powered_rate(r_MCC[MCC==3]),AER=get_AER(r_MCC[MCC==3])))
MCC_summary=rbind(MCC_summary,data.table(MCC=4, RR=get_RR(r_MCC[MCC==4]), RR_all=get_RR(r_MCC[MCC==4], include="all"), Power=get_powered_rate(r_MCC[MCC==4]),AER=get_AER(r_MCC[MCC==4])))
MCC_summary=rbind(MCC_summary,data.table(MCC=5, RR=get_RR(r_MCC[MCC==5]), RR_all=get_RR(r_MCC[MCC==5], include="all"), Power=get_powered_rate(r_MCC[MCC==5]),AER=get_AER(r_MCC[MCC==5])))
MCC_summary=rbind(MCC_summary,data.table(MCC=6, RR=get_RR(r_MCC[MCC==6]), RR_all=get_RR(r_MCC[MCC==6], include="all"), Power=get_powered_rate(r_MCC[MCC==6]),AER=get_AER(r_MCC[MCC==6])))
MCC_summary=rbind(MCC_summary,data.table(MCC=7, RR=get_RR(r_MCC[MCC==7]), RR_all=get_RR(r_MCC[MCC==7], include="all"), Power=get_powered_rate(r_MCC[MCC==7]),AER=get_AER(r_MCC[MCC==7])))
MCC_summary=rbind(MCC_summary,data.table(MCC=8, RR=get_RR(r_MCC[MCC==8]), RR_all=get_RR(r_MCC[MCC==8], include="all"), Power=get_powered_rate(r_MCC[MCC==8]),AER=get_AER(r_MCC[MCC==8])))
#MCC_summary$MCC=as.factor(MCC_summary$MCC)

RR_plot=ggplot(MCC_summary)+geom_point(aes(x=MCC, y=RR),stat="identity",size=2.8,color="red")+geom_line(aes(x=MCC, y=RR),linetype = "dashed",color="red")+
  geom_point(aes(x=MCC,y=RR_all),size=2.8,stat="identity",color="blue")+geom_line(aes(x=MCC,y=RR_all),linetype = "dashed")+
  scale_y_continuous(labels = scales::percent,n.breaks=10)+theme_classic()+
  scale_x_continuous(n.breaks=8)+
  labs(x="Minimum code count",y="Replication rate")
RR_plot

#geom_point(aes(x=MCC,y=RR_all),size=2.8,stat="identity",color="blue")+geom_line(aes(x=MCC,y=RR_all),linetype = "dashed")+

RR_plot_power=ggplot(MCC_summary)+geom_point(aes(x=MCC, y=RR),stat="identity",size=2.8,color="red")+geom_line(aes(x=MCC, y=RR),linetype = "dashed",color="red")+
  scale_y_continuous(labels = scales::percent,n.breaks=10)+theme_classic()+
  scale_x_continuous(n.breaks=8)+
  labs(x="Minimum code count",y=expression(RR[Power]*""))
RR_plot_power

RR_plot_all=ggplot(MCC_summary)+geom_point(aes(x=MCC, y=RR_all),stat="identity",size=2.8,color="blue")+geom_line(aes(x=MCC, y=RR_all),linetype = "dashed",color="blue")+
  scale_y_continuous(labels = scales::percent,n.breaks=10)+theme_classic()+
  scale_x_continuous(n.breaks=8)+
  labs(x="Minimum code count",y=expression(RR[All]*""))
RR_plot_all

RR_plot_power=ggplot(MCC_summary)+geom_point(aes(x=MCC, y=Power),stat="identity",size=2.8,color="green")+geom_line(aes(x=MCC, y=Power),linetype = "dashed",color="green")+
  scale_y_continuous(labels = scales::percent,n.breaks=10)+theme_classic()+
  scale_x_continuous(n.breaks=8)+
  labs(x="Minimum code count",y="Percent powered")
RR_plot_power



table(r_MCC$MCC)

foo=r_MCC[powered==1,.(mean=mean(rOR),se=sd(rOR)/sqrt(.N)),by="MCC"]
foo$dataset=as.factor(foo$dataset)

#r_MCC$MCC=as.factor(r_MCC$MCC)
OR_plot=ggplot(foo,aes(x=MCC,y=mean))+geom_bar(stat="identity",fill="grey50")+coord_cartesian(ylim=c(1.1,1.3))+
  geom_errorbar(aes(ymax=mean + se, ymin=mean-se,width = 0.4))+theme_classic()+
  scale_y_continuous( labels = scales::number_format(accuracy = 0.01))+labs(y="Mean odds ratio",x="Minimum code count")
OR_plot
RR_plot
ggarrange(RR_plot,OR_plot,labels=c("A","B"))
ggsave("MCC.png")


compare_annotated_results(r_MCC[MCC==1],r_MCC[MCC==2])

RR_plot=ggplot(MCC_summary, aes(x=rand, y=RR))+geom_point(stat="identity",size=2.8)+geom_line(linetype = "dashed")+ scale_x_continuous(labels = scales::percent,n.breaks=10)+theme_classic()+
  labs(x="Percent randomized",y=expression(RR[Power]*""))

#get_AER(r_MCC1)
