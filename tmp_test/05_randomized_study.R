library(data.table)
library(gaston)
library(pgrm)
library(glue)
library(ggplot2)
library(ggpubr)
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

r060=run_rand_assoc(rand_frac=0.60)
write.table(r040,file="tmp_test/OUT/r060.csv",row.names = F, col.names = T)

r070=run_rand_assoc(rand_frac=0.70)
write.table(r070,file="tmp_test/OUT/r070.csv",row.names = F, col.names = T)

r080=run_rand_assoc(rand_frac=0.80)
write.table(r080,file="tmp_test/OUT/r080.csv",row.names = F, col.names = T)

r090=run_rand_assoc(rand_frac=0.90)
write.table(r090,file="tmp_test/OUT/r090.csv",row.names = F, col.names = T)


rand_summary=data.table()


rand_summary=data.table(rand=0, RR=get_RR(r000), AER=get_AER(r000))
rand_summary=rbind(rand_summary,data.table(rand=10, RR=get_RR(r010), AER=get_AER(r010)))

rand_summary=rbind(rand_summary,data.table(rand=20, RR=get_RR(r020), AER=get_AER(r020)))
rand_summary=rbind(rand_summary,data.table(rand=30, RR=get_RR(r030), AER=get_AER(r030)))
rand_summary=rbind(rand_summary,data.table(rand=40, RR=get_RR(r040), AER=get_AER(r040)))
rand_summary=rbind(rand_summary,data.table(rand=50, RR=get_RR(r050), AER=get_AER(r050)))
rand_summary=rbind(rand_summary,data.table(rand=60, RR=get_RR(r060), AER=get_AER(r060)))
rand_summary=rbind(rand_summary,data.table(rand=70, RR=get_RR(r070), AER=get_AER(r070)))
rand_summary=rbind(rand_summary,data.table(rand=80, RR=get_RR(r080), AER=get_AER(r080)))
rand_summary=rbind(rand_summary,data.table(rand=90, RR=get_RR(r090), AER=get_AER(r090)))
rand_summary=rbind(rand_summary,data.table(rand=100, RR=get_RR(r100), AER=get_AER(r100)))
write.table(rand_summary,"results_rand.txt",row.names = F,col.names = T,sep="\t")


###

rand_summary=read.table("results_rand.txt",header=TRUE)
rand_summary$rand=  rand_summary$rand/100
RR_plot=ggplot(rand_summary, aes(x=rand, y=RR))+geom_point(stat="identity",size=2.8)+geom_line(linetype = "dashed")+ scale_x_continuous(labels = scales::percent,n.breaks=10)+theme_classic()+
  labs(x="Percent randomized",y=expression(RR[Power]*""))
RR_plot
r000$dataset="0%"
r010$dataset="10%"
r020$dataset="20%"
r030$dataset="30%"
r040$dataset="40%"
r050$dataset="50%"
r060$dataset="60%"
r070$dataset="70%"
r080$dataset="80%"
r090$dataset="90%"
r100$dataset="100%"

rand_all = rbind(r000,r010)
rand_all = rbind(rand_all,r020)
rand_all = rbind(rand_all,r030)
rand_all = rbind(rand_all,r040)
rand_all = rbind(rand_all,r050)
rand_all = rbind(rand_all,r060)
rand_all = rbind(rand_all,r070)
rand_all = rbind(rand_all,r080)
rand_all = rbind(rand_all,r090)
rand_all = rbind(rand_all,r100)

foo=rand_all[,.(mean=mean(rOR),se=sd(rOR)/sqrt(.N)),by="dataset"]
foo$dataset=as.factor(foo$dataset)
OR_plot=ggplot(foo,aes(x=reorder(dataset,-mean),y=mean))+geom_bar(stat="identity",fill="grey50")+coord_cartesian(ylim=c(1,1.1))+
  geom_errorbar(aes(ymax=mean + se, ymin=mean-se,width = 0.4))+theme_classic()+
  scale_y_continuous( labels = scales::number_format(accuracy = 0.01))+ labs(y="Mean odds ratio",x="Percent randomized")
OR_plot
RR_plot
ggarrange(RR_plot,OR_plot,labels=c("A","B"))
ggsave("Randomize.png")

get_RR(r000)
get_RR(r010)
compare_annotated_results(r000,r010)
compare_annotated_results(r010,r020)
compare_annotated_results(r020,r030)
compare_annotated_results(r030,r040)
compare_annotated_results(r040,r050)
compare_annotated_results(r050,r060)
compare_annotated_results(r060,r070)
compare_annotated_results(r070,r080)
compare_annotated_results(r090,r100)


# r100$rand= '100%'
# r075$rand= '75%'
# r050$rand= '50%'
# r025$rand= '25%'
# r000$rand= '0%'
#
#
#
# r075$AsubE= r075$rep-r075$Power
# r050$AsubE= r050$rep-r050$Power
# r025$AsubE= r025$rep-r025$Power
# r000$AsubE= (r000$rep)-(r000$Power)
# r100$AsubE= (r100$rep-1)/(r100$Power)
# mean(r100$AsubE)
# sum(r100$rep)/sum(r100$Power)
# get_RR(r100)
#
#
# r_all=rbind(r100,r075)
# r_all=rbind(r_all,r025)
# r_all=rbind(r_all,r050)
# r_all=rbind(r_all,r000)
#
# r_all$AER = r_all$rep/r_all$Power
#
# hist(r_all$AER)
#
# m        = rev(c(mean(r000$AsubE), mean(r075$AsubE),mean(r050$AsubE),mean(r025$AsubE),mean(r000$AsubE)))
# names(m) = rev(c("100%", "75%","50%","25%","0%"))
# se       = rev(c(sd(r100$AsubE)/sqrt(length(r100$AsubE)),
#              sd(r075$AsubE)/sqrt(length(r075$AsubE)),
#              sd(r050$AsubE)/sqrt(length(r050$AsubE)),
#              sd(r025$AsubE)/sqrt(length(r025$AsubE)),
#              sd(r000$AsubE)/sqrt(length(r000$AsubE))))
#
# bp= barplot(m, ylim=c(1,1.1), xpd=FALSE)
# box()
# arrows(x0=bp, y0=m-se, y1=m+se, code=3, angle=90)
#
# library(ggplot2)
# ggplot(r_all,aes(x=rand, y=rOR))+geom_boxplot()+coord_cartesian(ylim=c(.9,1.4))
#
#
# m        = rev(c(mean(r100$rOR), mean(r075$rOR),mean(r050$rOR),mean(r025$rOR),mean(r000$rOR)))
# names(m) = rev(c("100%", "75%","50%","25%","0%"))
# se       = rev(c(sd(r100$rOR)/sqrt(length(r100$rOR)),
#              sd(r075$rOR)/sqrt(length(r075$rOR)),
#              sd(r050$rOR)/sqrt(length(r050$rOR)),
#              sd(r025$rOR)/sqrt(length(r025$rOR)),
#              sd(r000$rOR)/sqrt(length(r000$rOR))
#              ))
#
# bp= barplot(m, ylim=c(1,1.1), xpd=FALSE)
# box()
# arrows(x0=bp, y0=m-se, y1=m+se, code=3, angle=90)
# bp
# m
#
# compare_annotated_results(r100,r075)
