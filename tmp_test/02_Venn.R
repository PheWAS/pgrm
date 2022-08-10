

table(benchmark_results$cohort)

d=benchmark_results
table(d$cohort)

studies_in_all=intersect(d[cohort=="BioVU_EUR" & powered==1]$assoc_ID,d[cohort=="UKBB" & powered==1]$assoc_ID)
studies_in_all=intersect(studies_in_all,d[cohort=="MGI" & powered==1]$assoc_ID)
length(studies_in_all) ## 393
d=d[assoc_ID %in% studies_in_all]

no_rep = intersect(d[cohort=="BioVU_EUR" & rep==0]$assoc_ID, d[cohort=="UKBB" & rep==0]$assoc_ID)
no_rep = intersect(no_rep, d[cohort=="MGI" & rep==0]$assoc_ID )
no_rep = intersect(studies_in_all,no_rep)
length(no_rep)

table(d[assoc_ID %in% no_rep & d$cohort == "UKBB"]$category_string)
library(formattable)

formattable(d[cohort=="UKBB" &  assoc_ID %in% no_rep,c('assoc_ID','rsID','phecode','phecode_string','category_string','Study_accession')],
            align =c("l","l","l","l","l"),style = "height:200px; overflow-y: scroll;overflow-x: scroll;")

table(d[cohort=="UKBB" &  assoc_ID %in% no_rep]$Study_accession)
table(d[cohort=="BioVU_EUR"]$category_string)

6+7+9+42+258+32+16


x <- list (
  BioVU =  d[cohort=="BioVU_EUR" & powered==1 & rep==1 ]$assoc_ID,
  MGI =  d[cohort=="MGI" & powered==1 & rep==1]$assoc_ID,
  UKBB =  d[cohort=="UKBB" & powered==1 & rep==1]$assoc_ID
  #,null =  d[d$null==1,]$assoc_ID
)
library(ggvenn)
ggvenn(
  x,
  fill_color = c("#EFC000FF", "#0073C2FF", "#CD534CFF", "#CD534CFF"),
  stroke_size = 0.8, set_name_size = 8, show_percentage=TRUE,digits=0,text_size=4.2
)


258 ## rep in all
258+42+32+7 ## rep in 2+

d[d$assoc_ID %in% studies_in_all &
    !assoc_ID %in% d[cohort=="BioVU_EUR" & powered==1 & rep==1 ]$assoc_ID &
    !assoc_ID %in% d[cohort=="MGI" & powered==1 & rep==1]$assoc_ID &
    !assoc_ID %in% d[cohort=="UKBB" & powered==1 & rep==1]$assoc_ID ]

#5+7+10+41+258+31+16
# 368
