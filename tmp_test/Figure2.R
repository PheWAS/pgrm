library(ggplot2)
library(ggpubr)

rand_summary=read.table("results_rand.txt",header=TRUE)
rand_summary$rand=  rand_summary$rand/100
fig_2a=ggplot(rand_summary, aes(x=rand, y=RR))+geom_point(stat="identity",size=2.8)+geom_line(linetype = "dashed")+ scale_x_continuous(labels = scales::percent,n.breaks=10)+theme_classic()+
  labs(x="Percent randomized",y=expression(RR[Power]*""))+scale_y_continuous(labels = scales::percent) ## ,n.breaks=5,limits=c(0,1)
fig_2a


## inpatient only

r_inpt_only=fread(file="anno_BioVU_EUR_INPT_only.csv",header=T,colClasses = list(character = 'phecode'))


compare_annotated_results(benchmark_results[cohort=="BioVU_EUR"],r_inpt_only)

fig2=data.frame(set="Inpatient only",odds_ratio=1/0.9895, U95=1/0.7444, L95=1/1.332)

r_no_exclude=fread(file="anno_BioVU_EUR_no_exclude.csv",header=T,colClasses = list(character = 'phecode'))

compare_annotated_results(r_no_exclude,benchmark_results[cohort=="BioVU_EUR"])

fig2=rbind(fig2,
           data.frame(set="No exclude range", odds_ratio=1/0.9294, U95=1/0.7453, L95=1/1.1591))
#fig2=rbind(fig2,data.frame(set="no exclude range", odds_ratio=1.07, L95=0.86, U95=1.34))

fig2

ggplot(fig2,aes(y=set, x=odds_ratio,xmin=L95, xmax=U95))+geom_pointrange()+geom_vline(xintercept=1, lty=2)+
  theme_classic()+labs(x="Odds ratio",y="")

fig_2b=ggplot(fig2,aes(y=set, x=odds_ratio))+geom_point(size=2)+geom_errorbar(aes(xmin=L95,xmax=U95,width=0.1))+geom_vline(xintercept=1, lty=2)+
  theme_classic()+labs(x="Odds ratio",y="")+scale_x_continuous(limits=c(0.5,1.5))
fig_2b
## MCC

MCC2=fread("~/Documents/GitHub/pgrm/tmp_test/OUT/results_BioVU_EUR_2.csv")
MCC2=MCC2[,c("assoc_ID")]
r_MCC=data.table()
for(i in 1:8){
  infile=glue("~/Documents/GitHub/pgrm/tmp_test/OUT/results_BioVU_EUR_{i}.csv")
  r=fread(infile)
  r=merge(MCC2,r,by="assoc_ID",all.x=T)
  r[is.na(r)]=0
  r$MCC= i
  r_MCC=rbind(r_MCC,r)
  #write.table(r,file=outfile,row.names = F, col.names = T)
}

MCC_summary=data.table(MCC=1, RR=get_RR(r_MCC[MCC==1]), assoc_tested=nrow(r_MCC[MCC==1]),RR_all=get_RR(r_MCC[MCC==1], include="all"),nPower=sum(r_MCC[MCC==1]$powered), Power=get_powered_rate(r_MCC[MCC==1]),AER=get_AER(r_MCC[MCC==1]))
MCC_summary=rbind(MCC_summary,data.table(MCC=2, assoc_tested=nrow(r_MCC[MCC==2]),RR=get_RR(r_MCC[MCC==2]), RR_all=get_RR(r_MCC[MCC==2], include="all"), nPower=sum(r_MCC[MCC==2]$powered), Power=get_powered_rate(r_MCC[MCC==2]), AER=get_AER(r_MCC[MCC==2])))
MCC_summary=rbind(MCC_summary,data.table(MCC=3, assoc_tested=nrow(r_MCC[MCC==3]),RR=get_RR(r_MCC[MCC==3]), RR_all=get_RR(r_MCC[MCC==3], include="all"), nPower=sum(r_MCC[MCC==3]$powered),Power=get_powered_rate(r_MCC[MCC==3]),AER=get_AER(r_MCC[MCC==3])))
MCC_summary=rbind(MCC_summary,data.table(MCC=4,assoc_tested=nrow(r_MCC[MCC==4]), RR=get_RR(r_MCC[MCC==4]), RR_all=get_RR(r_MCC[MCC==4], include="all"), nPower=sum(r_MCC[MCC==4]$powered),Power=get_powered_rate(r_MCC[MCC==4]),AER=get_AER(r_MCC[MCC==4])))
MCC_summary=rbind(MCC_summary,data.table(MCC=5, assoc_tested=nrow(r_MCC[MCC==5]),RR=get_RR(r_MCC[MCC==5]), RR_all=get_RR(r_MCC[MCC==5], include="all"), nPower=sum(r_MCC[MCC==5]$powered),Power=get_powered_rate(r_MCC[MCC==5]),AER=get_AER(r_MCC[MCC==5])))
MCC_summary=rbind(MCC_summary,data.table(MCC=6,assoc_tested=nrow(r_MCC[MCC==6]), RR=get_RR(r_MCC[MCC==6]), RR_all=get_RR(r_MCC[MCC==6], include="all"), nPower=sum(r_MCC[MCC==6]$powered),Power=get_powered_rate(r_MCC[MCC==6]),AER=get_AER(r_MCC[MCC==6])))
MCC_summary=rbind(MCC_summary,data.table(MCC=7, assoc_tested=nrow(r_MCC[MCC==7]),RR=get_RR(r_MCC[MCC==7]), RR_all=get_RR(r_MCC[MCC==7], include="all"), nPower=sum(r_MCC[MCC==7]$powered),Power=get_powered_rate(r_MCC[MCC==7]),AER=get_AER(r_MCC[MCC==7])))
MCC_summary=rbind(MCC_summary,data.table(MCC=8, assoc_tested=nrow(r_MCC[MCC==8]),RR=get_RR(r_MCC[MCC==8]), RR_all=get_RR(r_MCC[MCC==8], include="all"), nPower=sum(r_MCC[MCC==8]$powered),Power=get_powered_rate(r_MCC[MCC==8]),AER=get_AER(r_MCC[MCC==8])))
MCC_summary

r_MCC$powered=as.factor(r_MCC$powered)
r_MCC$MCC_fac=as.factor(r_MCC$MCC)
r_MCC$MCC_fac=relevel(r_MCC$MCC_fac,ref="2")

m_RRp=glm(rep~MCC_fac,data=r_MCC[powered==1],family="binomial")
m_RRa=glm(rep~MCC_fac,data=r_MCC,family="binomial")

m_power=glm(powered~MCC_fac,data=r_MCC,family="binomial")


p_RRa=predict(m_RRa,newdata=data.frame(MCC_fac=rep(unique(r_MCC$MCC_fac))),type="response",se.fit=T)
p_RRa=data.frame(p_RRa)
p_RRa$MCC=rownames(p_RRa)

plot_RR_ALL=ggplot(p_RRa,aes(x=MCC,y=fit))+geom_point(size=2.4,color="blue")+geom_errorbar(aes(x=MCC,ymin=fit-se.fit,ymax=fit+se.fit,width=0.2,fill="blue"))+
  theme_classic()+scale_y_continuous(labels = scales::percent,limits=c(0.30,.45))+labs(x="Minimum code count",y=expression(RR[All]*""))
plot_RR_ALL

p_RRp=predict(m_RRp,newdata=data.frame(MCC_fac=rep(unique(r_MCC$MCC_fac))),type="response",se.fit=T)
p_RRp=data.frame(p_RRp)
p_RRp$MCC=rownames(p_RRp)
p_RRp=merge(p_RRp,MCC_summary,by="MCC")



plot_RR_power=ggplot(p_RRp,aes(x=MCC,y=fit,label=nPower))+geom_point(size=2.4,color="red")+geom_errorbar(aes(x=MCC,ymin=fit-se.fit,ymax=fit+se.fit,width=0.2))+
  theme_classic()+scale_y_continuous( labels = scales::percent,limits=c(0.65,.85))+labs(x="Minimum code count",y=expression(RR[Power]*""))+geom_text(hjust=.5, vjust=-2)
plot_RR_power

# RR_plot_powered=ggplot(MCC_summary,aes(x=MCC, y=RR,label=nPower))+geom_point(,stat="identity",size=2.8,color="red")+geom_line(aes(x=MCC, y=RR),linetype = "dashed",color="red")+
#   scale_y_continuous(labels = scales::percent_format(accuracy = 5L),limits=c(0.65,.85))+theme_classic()+
#   scale_x_continuous(n.breaks=8)+geom_text(hjust=.4, vjust=-1)+
#   labs(x="Minimum code count",y=expression(RR[Power]*""))
# RR_plot_powered
#
# RR_plot_all=ggplot(MCC_summary)+geom_point(aes(x=MCC, y=RR_all),stat="identity",size=2.8,color="blue")+geom_line(aes(x=MCC, y=RR_all),linetype = "dashed",color="blue")+
#   scale_y_continuous(labels = scales::percent_format(accuracy = 5L),limits=c(0.35,.42))+theme_classic()+
#   scale_x_continuous(n.breaks=8)+
#   labs(x="Minimum code count",y=expression(RR[All]*""))
# RR_plot_all
#
# RR_plot_power=ggplot(MCC_summary)+geom_point(aes(x=MCC, y=Power),stat="identity",size=2.8,color="green")+geom_line(aes(x=MCC, y=Power),linetype = "dashed",color="green")+
#   scale_y_continuous(labels = scales::percent,n.breaks=10)+theme_classic()+
#   scale_x_continuous(n.breaks=8)+
#   labs(x="Minimum code count",y="Percent powered")
# RR_plot_power
fig_2c=plot_RR_ALL
fig_2d=plot_RR_power


ggarrange(fig_2a,fig_2b,fig_2c,fig_2d,labels=c("A","B","C","D"))
ggsave("Fig2_tmp.png")




TRP_C<-100/(100+650)
FPR_C<-200/(200+650)
C<-data.frame(TPR=TRP_C,FPR=FPR_C)

TRP_D<-120/(120+30)
FPR_D<-350/(350+500)
D<-data.frame(TPR=TRP_D,FPR=FPR_D)

lt_red='#f1b6da'
dk_red='#d01c8b'
lt_green='#b8e186'
dk_green='#4dac26'
d=data.frame(study="No Exclude Range",RR_ALL=lt_red,RR_power=lt_red,Power=lt_green)


d=data.frame(study="No Exclude Range",measure="RR_All",c='Nominally less')
d=rbind(d,data.frame(study="No Exclude Range",measure="RR_power",c='Nominally less'))
d=rbind(d,data.frame(study="No Exclude Range",measure="Powered",c='Nominally greater'))

d=rbind(d,data.frame(study="Inpatient only",measure="RR_All",c='Significantly less'))
d=rbind(d,data.frame(study="Inpatient only",measure="RR_power",c='Nominally less'))
d=rbind(d,data.frame(study="Inpatient only",measure="Powered",c='Significantly greater'))

d=rbind(d,data.frame(study="MCC=1",measure="RR_All",c='Nominally less'))
d=rbind(d,data.frame(study="MCC=1",measure="RR_power",c='Significantly less'))
d=rbind(d,data.frame(study="MCC=1",measure="Powered",c='Significantly greater'))

d=rbind(d,data.frame(study="MCC=3",measure="RR_All",c='Nominally less'))
d=rbind(d,data.frame(study="MCC=3",measure="RR_power",c='Nominally greater'))
d=rbind(d,data.frame(study="MCC=3",measure="Powered",c='Significantly less'))

d=rbind(d,data.frame(study="MCC=4",measure="RR_All",c='Significantly less'))
d=rbind(d,data.frame(study="MCC=4",measure="RR_power",c='Nominally greater'))
d=rbind(d,data.frame(study="MCC=4",measure="Powered",c='Significantly less'))

d=rbind(d,data.frame(study="MCC=5",measure="RR_All",c='Significantly less'))
d=rbind(d,data.frame(study="MCC=5",measure="RR_power",c='Significantly greater'))
d=rbind(d,data.frame(study="MCC=5",measure="Powered",c='Significantly less'))

d$measure=as.factor(d$measure)
d$measure <- factor(d$measure, levels=c("RR_All","Powered","RR_power"))

d$study <- factor(d$study, levels=rev(c("No Exclude Range","Inpatient only","MCC=1","MCC=3","MCC=4","MCC=5")))
ggplot(d,aes(x=measure,y=study, color= c) )+geom_point(size=20)+
  scale_color_manual(breaks = c("Nominally less", "Significantly less", "Nominally greater","Significantly greater"),values=c(lt_red,dk_red,lt_green,dk_green))+
  theme_classic()


summary(m_RRa)
summary(m_RRp)
summary(m_power)
