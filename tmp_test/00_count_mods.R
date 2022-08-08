library('data.table')
library(pgrm)
library(effects)
rawDir = 'tmp_test'

PGRM_mod = fread(
  file.path(rawDir, 'PGRM_mods_info.csv'),
  colClasses = list(character = 'phecode'),
  na.strings = '-')

PGRM_mod$study_date=as.Date(PGRM_mod$study_date)
min(PGRM_mod$study_date)
PGRM_mod$study_day=as.numeric(difftime(PGRM_mod$study_date,min(PGRM_mod$study_date),units = "days"),units = "days")
hist(PGRM_mod$study_day)
PGRM_mod$mod_pheno=as.factor(PGRM_mod$mod_pheno)
PGRM_mod$backgrou=as.factor(PGRM_mod$backgrou)
PGRM_mod$alt_model=as.factor(PGRM_mod$alt_model)

PGRM_mod[,.(mod_pheno,backgrou,alt_model),by=pubmedid]

m=glm(study_day~mod_pheno+backgrou+alt_model, data=PGRM_mod[,.(mod_pheno,backgrou,alt_model,study_day),by=pubmedid])
summary(m)
plot(allEffects(m))
table(PGRM_mod$mod_pheno)
table(PGRM_mod$backgrou)
table(PGRM_mod$alt_model)

table(PGRM_mod$modified_pheno)
prop.table(table(PGRM_mod$modified_pheno))
