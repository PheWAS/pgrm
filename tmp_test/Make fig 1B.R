
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(wesanderson)
d=get_PGRM(unique=F)
table(d$ancestry)
uniq=get_PGRM(unique=T)
nrow(d)
nrow(uniq)

counts=uniq[,.(n=.N), by=list(phecode,phecode_string)]
phecodes_include = counts[n>=50,]$phecode
length(phecodes_include)

#foo=d[,.(n=.N,ancestry=paste(unique(ancestry), collapse=', ')),by=list(phecode, phecode_string,category_string )]
#foo
#d=d[,.(ancestry=paste(unique(ancestry), collapse=', ')),by=list(SNP,phecode, phecode_string,category_string )]
#d[grep(",", ancestry)]$ancestry = 'Multiple'
#table(d$ancestry)

d[phecode=="731.1"]$phecode_string="Osteitis deformans"
d[phecode=="362.29"]$phecode_string="Age related macular degeneration"
d[phecode=="428"]$phecode_string="Congestive heart failure"
d[phecode=="523.3"]$phecode_string="Periodontitis"

d[phecode=="728.71"]$phecode_string="Dupuytren's disease"
d[phecode=="443.9"]$phecode_string="Peripheral vascular disease"
d[phecode=="189.11"]$phecode_string="Renal cancer"
d[phecode=="189.21"]$phecode_string="Bladder cancer"
d[phecode=="185"]$phecode_string="Prostate cancer"
d[phecode=="174.11"]$phecode_string="Breast cancer"
d[phecode=="159.3"]$phecode_string="Gallbladder cancer"
d[phecode=="184.11"]$phecode_string="Ovarian cancer"
d[phecode=="153.3"]$phecode_string="Rectal cancer"
d[phecode=="172.1"]$phecode_string="Melanoma"

d[phecode=="571.51"]$phecode_string="Cirrhosis w/out mention of alcohol"
d[phecode=="585.3"]$phecode_string="Chronic renal failure"
d[phecode=="165.1"]$phecode_string="Lung cancer"
d[phecode=="xxxx"]$phecode_string="xxxx"


d[phecode=="189.11"]

length(phecodes_include) ## 28
d$ancestry = factor(d$ancestry, levels=rev(c('EUR','EAS','AFR','AMR','SAS')))
table(d$ancestry)
p<-ggplot(data=d[phecode %in% phecodes_include,], aes(x=reorder(phecode_string, phecode_string, function(x)-length(x)*-1),fill=ancestry)) +
  geom_bar(color="grey40")+scale_y_continuous(expand = c(0, 0))+scale_fill_brewer(palette="Set2", direction=-1)
p=p+ coord_flip(expand=T)+theme_classic()+xlab("")+ylab("Number of SNPs")
  #scale_fill_brewer(palette="Set2", direction=-1)
#+scale_fill_manual(values=rev(c( "#377eb8","#e41a1c", "#4daf4a","#984ea3","#ff7f00")))
p
table(d$category_string)
ggsave(p,file="figure_s1.png") ## Saving 11.2 x 6.82 in image

d[category_string=="metabolic/heme"]$category_string="metab."
d[category_string=="respiratory"]$category_string="resp."
d[category_string=="infectious diseases"]$category_string="infectious"
d[category_string=="musculoskeletal"]$category_string="musculoskel."
d[category_string=="genitourinary"]$category_string="genitourinary"
p<-ggplot(data=d[!phecode %in% phecodes_include,], aes(x=reorder(phecode_string, phecode_string, function(x)-length(x)*-1),fill=ancestry)) +
  geom_bar(color="grey40")+scale_y_continuous(expand = c(0, 0))+scale_fill_brewer(palette="Set2", direction=-1)
p=p + facet_grid(category_string~.,scales="free_y", space = "free")+coord_flip()+theme_classic()+xlab("")+ylab("Number of SNPs")
ggsave(p,file="figure_s2.png") ## Saving 11.6 x 14.3 in image

nrow(d[!phecode %in% phecodes_include,])##1307
nrow(d[!phecode %in% phecodes_include,])##1307
length(counts[n<50,]$phecode)
