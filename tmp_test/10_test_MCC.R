##foobar
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
MCC_summary=data.table(MCC=1, RR=get_RR(r_MCC[MCC==1]), assoc_tested=.N,RR_all=get_RR(r_MCC[MCC==1], include="all"),nPower=sum(r_MCC[MCC==1]$powered), Power=get_powered_rate(r_MCC[MCC==1]),AER=get_AER(r_MCC[MCC==1]))
MCC_summary=rbind(MCC_summary,data.table(MCC=2, assoc_tested=.N,RR=get_RR(r_MCC[MCC==2]), RR_all=get_RR(r_MCC[MCC==2], include="all"), nPower=sum(r_MCC[MCC==2]$powered), Power=get_powered_rate(r_MCC[MCC==2]), AER=get_AER(r_MCC[MCC==2])))
MCC_summary=rbind(MCC_summary,data.table(MCC=3, assoc_tested=.N,RR=get_RR(r_MCC[MCC==3]), RR_all=get_RR(r_MCC[MCC==3], include="all"), nPower=sum(r_MCC[MCC==3]$powered),Power=get_powered_rate(r_MCC[MCC==3]),AER=get_AER(r_MCC[MCC==3])))
MCC_summary=rbind(MCC_summary,data.table(MCC=4,assoc_tested=.N, RR=get_RR(r_MCC[MCC==4]), RR_all=get_RR(r_MCC[MCC==4], include="all"), nPower=sum(r_MCC[MCC==4]$powered),Power=get_powered_rate(r_MCC[MCC==4]),AER=get_AER(r_MCC[MCC==4])))
MCC_summary=rbind(MCC_summary,data.table(MCC=5, assoc_tested=.N,RR=get_RR(r_MCC[MCC==5]), RR_all=get_RR(r_MCC[MCC==5], include="all"), nPower=sum(r_MCC[MCC==5]$powered),Power=get_powered_rate(r_MCC[MCC==5]),AER=get_AER(r_MCC[MCC==5])))
MCC_summary=rbind(MCC_summary,data.table(MCC=6,assoc_tested=.N, RR=get_RR(r_MCC[MCC==6]), RR_all=get_RR(r_MCC[MCC==6], include="all"), nPower=sum(r_MCC[MCC==6]$powered),Power=get_powered_rate(r_MCC[MCC==6]),AER=get_AER(r_MCC[MCC==6])))
MCC_summary=rbind(MCC_summary,data.table(MCC=7, assoc_tested=.N,RR=get_RR(r_MCC[MCC==7]), RR_all=get_RR(r_MCC[MCC==7], include="all"), nPower=sum(r_MCC[MCC==7]$powered),Power=get_powered_rate(r_MCC[MCC==7]),AER=get_AER(r_MCC[MCC==7])))
MCC_summary=rbind(MCC_summary,data.table(MCC=8, assoc_tested=.N,RR=get_RR(r_MCC[MCC==8]), RR_all=get_RR(r_MCC[MCC==8], include="all"), nPower=sum(r_MCC[MCC==8]$powered),Power=get_powered_rate(r_MCC[MCC==8]),AER=get_AER(r_MCC[MCC==8])))
MCC_summary
#MCC_summary$MCC=as.factor(MCC_summary$MCC)

compare_annotated_results(r_MCC[MCC==1],r_MCC[MCC==2])
compare_annotated_results(r_MCC[MCC==3],r_MCC[MCC==2])

compare_annotated_results(r_MCC[MCC==1],r_MCC[MCC==2],include_all = F)
compare_annotated_results(r_MCC[MCC==2],r_MCC[MCC==3],include_all = F)
compare_annotated_results(r_MCC[MCC==3],r_MCC[MCC==4],include_all = F)
compare_annotated_results(r_MCC[MCC==4],r_MCC[MCC==5],include_all = F)
compare_annotated_results(r_MCC[MCC==5],r_MCC[MCC==6],include_all = F)
compare_annotated_results(r_MCC[MCC==6],r_MCC[MCC==7],include_all = F)
compare_annotated_results(r_MCC[MCC==7],r_MCC[MCC==8],include_all = F)

compare_annotated_results(r_MCC[MCC==3],r_MCC[MCC==6],include_all = F)

r_MCC$powered=as.factor(r_MCC$powered)
r_MCC$MCC_fac=as.factor(r_MCC$MCC)
r_MCC$MCC_fac=relevel(r_MCC$MCC_fac,ref="2")

m_RRp=glm(rep~MCC_fac,data=r_MCC[powered==1],family="binomial")
summary(m_RRp)
plot(allEffects(m_RRp))

m_RRa=glm(rep~MCC_fac,data=r_MCC,family="binomial")
summary(m_RRa)
plot(allEffects(m_RRa))

m_pow=glm(powered~m_pow,data=r_MCC,family="binomial")
summary(m_pow)
plot(allEffects(m_pow))

p_RRa=predict(m_RRa,newdata=data.frame(MCC_fac=rep(unique(r_MCC$MCC_fac))),type="response",se.fit=T)
p_RRa=data.frame(p_RRa)
p_RRa$MCC=rownames(p_RRa)

plot_RR_ALL=ggplot(p_RRa,aes(x=MCC,y=fit))+geom_point()+geom_errorbar(aes(x=MCC,ymin=fit-se.fit,ymax=fit+se.fit,width=0.4))+
  theme_classic()+scale_y_continuous( limits=c(.30,.45),labels = scales::percent)+labs(x="Minimum code count",y=expression(RR[All]*""))
plot_RR_ALL

p_RRp=predict(m_RRp,newdata=data.frame(MCC_fac=rep(unique(r_MCC$MCC_fac))),type="response",se.fit=T)
p_RRp=data.frame(p_RRp)
p_RRp$MCC=rownames(p_RRp)

plot_RR_ALL=ggplot(p_RRp,aes(x=MCC,y=fit))+geom_point()+geom_errorbar(aes(x=MCC,ymin=fit-se.fit,ymax=fit+se.fit,width=0.4))+
  theme_classic()+scale_y_continuous( limits=c(.30,.45),labels = scales::percent)+labs(x="Minimum code count",y=expression(RR[All]*""))
plot_RR_ALL


cat_strings= unique(r_MCC$category_string)
i=0

i=i+1
print(cat_strings[i])

compare_annotated_results(r_MCC[MCC==3 & category_string==cat_strings[i]],
                          r_MCC[MCC==4 & category_string==cat_strings[i]])
r_MCC$MCC_fac=as.factor(r_MCC$MCC)
r_MCC$MCC_fac=relevel(r_MCC$MCC_fac,ref="1")
m=glm(rep~MCC_fac:category_string,data=r_MCC[powered==1],family="binomial")
summary(m)
p=predict(m,newdata=data.frame(MCC_fac=rep(unique(r_MCC$MCC_fac),length(unique(r_MCC$category_string))),
                               category_string=rep(unique(r_MCC$category_string),8)),type="response",se.fit=T)
p=data.frame()
cat_strings= unique(r_MCC$category_string)
for(i in 1:length(cat_strings)){
  pred=predict(m,newdata=data.frame(MCC_fac=unique(r_MCC$MCC_fac),category_string=rep(cat_strings[i],8)),type="response",se.fit=T)
  pred=data.frame(pred)
  pred$MCC=rownames(pred)
  pred$cat = cat_strings[i]
  head(pred)
  p=rbind(p,pred)
}
p$MCC=as.numeric(p$MCC)
library( ggrepel)
ggplot(p,aes(x=MCC,y=fit,color=cat))+geom_point()+geom_line()+ geom_text_repel(
  aes(label = cat), data = p,
  fontface ="plain", color = "black", size = 3
)

r_MCC$MCC_fac=relevel(r_MCC$MCC_fac,ref="2")
m=glm(rep~MCC_fac,data=r_MCC,family="binomial")
summary(m)
p=predict(m,newdata=data.frame(MCC_fac=unique(r_MCC$MCC_fac),category_string=rep('neoplasms',8)),type="response",se.fit=T)
p
foo=data.frame(p)
foo$MCC=rownames(foo)
plot_RR_ALL=ggplot(foo,aes(x=MCC,y=fit))+geom_point()+geom_errorbar(aes(x=MCC,ymin=fit-se.fit,ymax=fit+se.fit,width=0.4))+
  theme_classic()+scale_y_continuous( limits=c(.30,.45),labels = scales::percent)+labs(x="Minimum code count",y=expression(RR[All]*""))
plot_RR_ALL
#p=predict(m,newdata=data.frame(MCC_fac=unique(r_MCC$MCC_fac),category_string=rep('sense organs',8)),type="response",se.fit=T)

m=glm(rep~MCC_fac,data=r_MCC[powered==1],family="binomial")
summary(m)
p=predict(m,newdata=data.frame(MCC_fac=unique(r_MCC$MCC_fac),Power=rep(.8,8)),type="response",se.fit=T)
p
foo2=data.frame(p)
foo2$MCC=rownames(foo2)
plot_RR_power=ggplot(foo2,aes(x=MCC,y=fit))+geom_point()+geom_errorbar(aes(x=MCC,ymin=fit-se.fit,ymax=fit+se.fit,width=0.4))+
  theme_classic()+scale_y_continuous( labels = scales::percent)+labs(x="Minimum code count",y=expression(RR[Powered]*""))
plot_RR_power
ggarrange(plot_RR_ALL,plot_RR_power,labels="AUTO")





predict(m,newdata=data.frame(MCC_fac=unique(r_MCC$MCC_fac),category_string="sense organs"))

library(effects)
plot(allEffects(m))




#scale_y_continuous(labels = scales::percent,n.breaks=10,name = expression(RR[Power]*""),sec.axis = sec_axis( trans=~.*10000, name="Second Axis"))+
RR_plot=ggplot()+geom_bar(data=MCC_summary,aes(x=MCC,y=nPower),stat="identity")+
  geom_point(data=MCC_summary,aes(x=MCC, y=RR*100),stat="identity",size=2.8,color="red")+
  scale_y_continuous(name = "Powered associations",sec.axis = sec_axis(~ .*100,name= expression(RR[Power]*"")))+
  theme_classic()+
  scale_x_continuous(n.breaks=8)+
  labs(x="Minimum code count")
RR_plot
#geom_point(aes(x=MCC,y=RR_all),size=2.8,stat="identity",color="blue")+geom_line(aes(x=MCC,y=RR_all),linetype = "dashed")+

+geom_line(aes(x=MCC, y=RR),linetype = "dashed",color="red")+

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
