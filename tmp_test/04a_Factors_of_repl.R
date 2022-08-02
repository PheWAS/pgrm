library(pgrm)
library(splines)

d=benchmark_results

d$ancestry = as.factor(d$ancestry)
d$ancestry = relevel(d$ancestry, ref="EUR")

d$category_string = as.factor(d$category_string)
d$category_string = relevel(d$category_string, ref="neoplasms")

d$cohort = as.factor(d$cohort)
d$cohort = relevel(d$cohort, ref="MGI")

nrow(d) ## 10038

d$cat_L95=log(d$cat_L95)
d$first_pub_year=as.numeric(format(d$first_pub_date,"%Y"))


m=glm(rep~cat_OR+cat_LOG10_P+AF+cases+category_string+ancestry+first_pub_date, data=d,family="binomial")
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
