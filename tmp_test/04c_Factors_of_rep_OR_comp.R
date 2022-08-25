library(pgrm)
library(glue)
#library(splines)
#library(sjPlot)
#library(DescTools)
#library(ggplot2)
#library(jtools)

d=benchmark_results[rep==1]

d$ancestry = as.factor(d$ancestry)
d$ancestry = relevel(d$ancestry, ref="EUR")

nrow(d)

plot(d$rOR, d$cat_OR)

ggplot(d,aes(x=cat_OR,y=rOR))+geom_point()+theme_bw()+geom_smooth(method=lm)

ggplot(d[cat_OR<6],aes(x=cat_OR,y=rOR,color=cohort))+geom_point()+theme_bw()+geom_smooth(method=lm)+coord_cartesian(xlim=c(1,6),ylim=c(1,6))

m=lm(cat_OR~rOR,data=d)
summary(m)$adj.r.squared
R2=0.3437713


label = c()

OR_lim=4
ggplot(d,aes(x=cat_OR,y=rOR))+geom_point(aes(fill=cohort,alpha=0.8),pch=21,size=4)+theme_bw()+geom_smooth(method=lm,color="grey40")+coord_cartesian(xlim=c(1,OR_lim),ylim=c(1,OR_lim),expand=F)+
  geom_text(data = data.frame(x = 1, y = 0, label = "test"), aes(x = OR_lim-.4, y = OR_lim-.2, label = "R2 = 0.344"),parse=F)+xlab("Catalog Odds Ratio")+ylab("Test cohort Odds Ratio")

  geom_label(label="Look at this!", x=2,y=2)

  geom_text(x = 2, y = 6, label = label,parse = TRUE)


ggplot(d[cat_L95<4 & rL95 < 4],aes(x=cat_L95,y=rL95,color=cohort))+geom_point()+theme_bw()+geom_smooth(method=lm)+coord_cartesian(xlim=c(1,4),ylim=c(1,4))


ggplot(d[cat_L95<OR_lim & rL95 < OR_lim],aes(x=cat_L95,y=rL95,color=cohort))+geom_point()+theme_bw()+geom_smooth(method=lm)+coord_cartesian(xlim=c(1,OR_lim),ylim=c(1,OR_lim))

table(d$category_string)

circ=d[category_string=="musculoskeletal" & cohort %in% c("BioVU_EUR","MGI","UKBB")]

circ$rank=rank(circ$assoc_ID, ties.method="min")

ggplot(circ,aes(x=rOR, y=rank,color=cohort,group=phecode_string))+geom_point()+geom_errorbar(aes(xmin=rL95, xmax=rU95), width=.1)


(r1 <- rank(x1 <- c(3, 1, 4, 15, 92)))
x2 <- c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5)
names(x2) <- letters[1:11]
(r2 <- rank(x2)) # ties are averaged

tMeth <- eval(formals(rank)$ties.method)
rx2 <- sapply(tMeth, function(M) rank(x2, ties.method=M))
cbind(x2, rx2)
## ties.method's does not matter w/o ties:
x <- sample(47)
rx <- sapply(tMeth, function(MM) rank(x, ties.method=MM))
stopifnot(all(rx[,1] == rx))
