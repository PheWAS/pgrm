library(pgrm)
library(glue)
#library(splines)
#library(sjPlot)
#library(DescTools)
#library(ggplot2)
#library(jtools)

d=benchmark_results

d$ancestry = as.factor(d$ancestry)
d$ancestry = relevel(d$ancestry, ref="EUR")

foo=benchmark_results[,.(.N,rep=sum(rep)),by="phecode_string"]
nrow(foo[rep>0])

foo=d[powered==1,.(rr=sum(rep)/.N,sum(powered)),by="category_string"]
foo[order(rr)]
get_RR(d)
d$category_string = as.factor(d$category_string)
d$category_string = relevel(d$category_string, ref="neurological")

d$cohort = as.factor(d$cohort)
d$cohort = relevel(d$cohort, ref="BioVU_EUR")

d$log_cat_OR=log(d$cat_OR)
nrow(d) ## 10038
#d$rep = as.factor(d$rep)

t=t.test(d$cat_OR, d$rOR, paired=T)
t ## (mean of the differences 0.08316576 (0.07720165-0.08912986))
t$p.value

nrow(d[rep==1])
t=t.test(d[rep==1]$cat_OR, d[rep==1]$rOR, paired=T)
t
t$p.value
nrow(d)
d=merge(d,PGRM_pubinfo[,c('assoc_ID','pub_date')],by="assoc_ID")

m=glm(rep~log(cat_OR)+cat_LOG10_P+RAF+cases+controls+pub_date+pub_count+category_string+ancestry, data=d,family="binomial")
summary(m)
Cstat(m)
plot_model(m)
#tab_model(m)
df_m=data.frame(summary(m)$coef)
df_m
df_m$var = rownames(df_m)
df_m
CI=confint.default(m)
CI_df=data.frame(CI)
CI_df$var = rownames(CI)

glm_res = merge(df_m, CI_df, by="var")
glm_res
write.table(glm_res,file="factors_of_rep_ALL.txt",sep="\t",row.names = F,col.names = T)

cats= unique(d$category_string)
for(i in 1:length(cats)){
  print(glue('{i} {cats[i]}'))
  d$in_cat = 0
  d[d$category_string== cats[i]]$in_cat = 1
  m=glm(rep~log(cat_OR)+cat_LOG10_P+RAF+cases+controls+pub_date+pub_count+in_cat,data=d, family="binomial")
  summary(m)
}




#d$cat_L95=log(d$cat_L95)
d$first_pub_year=as.numeric(format(d$first_pub_date,"%Y"))

m=glm(rep~first_pub_year, data=d,family="binomial")
summary(m)
plot(allEffects(m))
OR=exp(-0.127056)
predict(m,newdata=data.frame(first_pub_year=2005),type="response")
p2005=0.7474751/(1-0.7474751)
predict(m,newdata=data.frame(first_pub_year=2006),type="response")
p2006=0.7227483/(1-0.7227483)

predict(m,newdata=data.frame(first_pub_year=2020),type="response")
predict(m,newdata=data.frame(first_pub_year=2010))

m=glm(rep~log(cat_OR)+cat_LOG10_P+cases+RAF, data=d,family="binomial")
with(summary(m), 1 - deviance/null.deviance)
summary(m)$coef
summary(m)

plot(allEffects(m))
anova(m)

m=glm(rep~ancestry, data=d,family="binomial")
summary(m)

m=glm(rep~first_pub_date, data=d,family="binomial")
summary(m)
year_OR = exp((-3.570e-04*365.25))
year_OR
predict(m,newdata=data.frame(first_pub_date=as.Date('2007-01-01')),type="response")
d[year(first_pub_date)==2016]

m=glm(rep~log_cat_OR+cat_LOG10_P+RAF+cases+first_pub_date, data=d,family="binomial")
summary(m)
summary(m)$coef
plot(allEffects(m))
with(summary(m), 1 - deviance/null.deviance)
Cstat(m)
odds.ratio(m)
#0.7833833
yr_OR = exp(-0.0003422429*365.25)

ndata=with(d, data.table(first_pub_date = seq(min(first_pub_date), max(first_pub_date),length = 1000),
           log_cat_OR=rep(mean(d$log_cat_OR),times=1000),
           cat_LOG10_P=rep(mean(d$cat_LOG10_P),times=1000),
           RAF=rep(mean(d$RAF),times=1000),
           cases=rep(mean(d$cases),times=1000)))
ndata$fit = predict(m, newdata = ndata, type = 'response')
head(ndata)
head(ndata[year(first_pub_date)==2007])
plt <- ggplot(ndata, aes(x = first_pub_date, y = fit)) +  geom_line()
plt
plot(allEffects(m))
#  ndata <- bind_cols(ndata, setNames(as_tibble(predict(mod, ndata, se.fit = TRUE)[1:2]), c('fit_link','se_link')))
## create the interval and backtransform
#ndata <- mutate(ndata, fit_resp = ilink(fit_link), right_upr = ilink(fit_link + (2  se_link)), right_lwr = ilink(fit_link - (2  se_link))) ## show ndata</code>

effect_plot(m)


m=glm(rep~log(cat_OR)+cat_LOG10_P+RAF+cases+cohort, data=d,family="binomial")
summary(m)
plot_model(m)
tab_model(m)

foo=d[,.(.N,rep=sum(rep)/.N),by="category_string"]
get_RR(d,include = "all")
foo[order(rep)]

m=glm(rep~first_pub_date, data=d,family="binomial")
summary(m)
plot(allEffects(m))
plot_model(m)




m=glm(rep~cat_OR+cat_LOG10_P+RAF+cases+cohort, data=d,family="binomial")
summary(m)
plot_model(m)

m=glm(rep~log(cat_OR)*cat_LOG10_P+RAF+cases+category_string+first_pub_date, data=d,family="binomial")
summary(m)
plot_model(m)
#tab_model(m)
df_m=data.frame(summary(m)$coef)
df_m
df_m$var = rownames(df_m)
df_m
CI=confint.default(m)
CI_df=data.frame(CI)
CI_df$var = rownames(CI)

glm_res = merge(df_m, CI_df, by="var")
glm_res
write.table(glm_res,file="factors_of_rep.txt",sep="\t",row.names = F,col.names = T)




t=table(d[rep==1]$cohort, d[rep==1]$CI_overlap)
t
write.table(t,file="overlap.txt",sep="\t")

plot(d[rep==1]$cat_OR,d[rep==1]$rOR )
cor(d[rep==1]$cat_OR,d[rep==1]$rOR )
t=t.test(d[rep==1]$cat_OR,d[rep==1]$rOR,paired=T )
t$p.value

plot(d[rep==1]$cat_L95,d[rep==1]$rL95 )
cor(d[rep==1]$cat_L95,d[rep==1]$rL95 )
t=t.test(d[rep==1]$cat_L95,d[rep==1]$rL95,paired=T )
t$p.value
t
nrow(d[rep==1])
nrow(d[rep==1 & cat_OR > rOR ])
nrow(d[rep==1 & cat_OR < rOR ])
