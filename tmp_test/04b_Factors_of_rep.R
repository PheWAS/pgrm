library(pgrm)
library(splines)
library(sjPlot)
library(ROCR)
library(caret)

d=benchmark_results

d$ancestry = as.factor(d$ancestry)
d$ancestry = relevel(d$ancestry, ref="EUR")

d$category_string = as.factor(d$category_string)
d$category_string = relevel(d$category_string, ref="neoplasms")

d$cohort = as.factor(d$cohort)
d$cohort = relevel(d$cohort, ref="MGI")

nrow(d) ## 10038
head(d)
#d$cat_L95=log(d$cat_L95)
d$first_pub_year=as.numeric(format(d$first_pub_date,"%Y"))
head(d)
table(d$cohort)
d$cc_ratio = d$cases/d$controls

trctrl <- trainControl(method = "cv", number = 10, savePredictions=TRUE)
nb_fit <- train(factor(rep) ~cat_OR+cat_LOG10_P+cases+cc_ratio+first_pub_date+category_string, data = d, method = "rf", trControl=trctrl, tuneLength = 0)
nb_fit <- train(factor(rep) ~cat_OR, data = d, method = "rf", trControl=trctrl, tuneLength = 0)
summary(nb_fit)

trctrl <- trainControl(method = "cv", number = 10)
nb_fit <- train(factor(rep) ~cat_OR+cat_LOG10_P+cases+controls+category_string, data = d[cohort=="UKBB"], method = "glm", family="binomial",trControl=trctrl, tuneLength = 0)

nb_fit <- train(P ~cat_OR+cat_LOG10_P+cases+controls+category_string, data = d, method = "glm",trControl=trctrl, tuneLength = 0)
summary(nb_fit)
nb_fit

BBJ=d[cohort=="BBJ"]

pred = predict(nb_fit,data=BBJ,type="prob")
nrow(pred)
nrow(MGI)
MGI$pred = pred$`1`
head(MGI)


pred <- nb_fit$pred
pred$equal <- ifelse(pred$pred == pred$obs, 1,0)

eachfold <- pred %>%
  group_by(Resample) %>%
  summarise_at(vars(equal),
               list(Accuracy = mean))
eachfold

ggplot(data=eachfold, aes(x=Resample, y=Accuracy, group=1)) +
  geom_boxplot(color="maroon") +
  geom_point() +
  theme_minimal()

m=glm(rep~cat_OR+cat_LOG10_P+AF+cases+controls+category_string, data=d[cohort %in% c("BioVU_EUR")],family="binomial")
summary(m)


BIOVU_EUR=d[cohort %in% c("BioVU_EUR")]
BIOVU_EUR$prob <- predict(m,newdata=BIOVU,type="response")
table(BIOVU_EUR[prob>.8]$rep)
prop.table(table(BIOVU_EUR[prob>.8]$rep)) ## 0.7691057

BIOVU_AFR=d[cohort %in% c("BioVU_AFR")]
BIOVU_AFR$prob <- predict(m,newdata=BIOVU_AFR,type="response")
table(BIOVU_AFR[prob>.8]$rep)
prop.table(table(BIOVU_AFR[prob>.8]$rep)) ## 1.00

MGI=d[cohort %in% c("MGI")]
MGI$prob <- predict(m,newdata=MGI,type="response")
table(MGI[prob>.5]$rep)
prop.table(table(MGI[prob>.8]$rep)) ## 0.7381657

BBJ=d[cohort %in% c("BBJ")]
BBJ$prob <- predict(m,newdata=BBJ,type="response")
table(BBJ[prob>.8]$rep)
prop.table(table(BBJ[prob>.8]$rep)) ## 0.7241379

UKBB=d[cohort %in% c("UKBB")]
UKBB$prob <- predict(m,newdata=UKBB,type="response")
table(UKBB[prob>.8]$rep)
prop.table(table(UKBB[prob>.8]$rep)) ## 0.93761468
plot(UKBB$prob,UKBB$Power)


ctrl <- trainControl(method = "repeatedcv", number = 10, savePredictions = TRUE)
mod_fit <- train(Class ~ Age + ForeignWorker + Property.RealEstate + Housing.Own +
                   CreditHistory.Critical,  data=GermanCredit, method="glm", family="binomial",
                 trControl = ctrl, tuneLength = 5)
pred = predict(mod_fit, newdata=testing)
confusionMatrix(data=pred, testing$Class)

ctrl <- trainControl(method = "repeatedcv", number = 10, savePredictions = TRUE)



pred <- prediction(prob, d[cohort %in% c("BioVU_EUR")]$rep)

perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

auc <- performance(pred, measure = "auc")
auc <- auc@y.values[[1]]
auc

#calc_class_err

plot_model(m)
tab_model(m)
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
