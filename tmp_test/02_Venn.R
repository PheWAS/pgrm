

table(benchmark_results$cohort)

d=benchmark_results
table(d$cohort)

studies_in_all=intersect(d[cohort=="BioVU_EUR" & powered==1]$assoc_ID,d[cohort=="UKBB" & powered==1]$assoc_ID)
studies_in_all=intersect(studies_in_all,d[cohort=="MGI" & powered==1]$assoc_ID)
length(studies_in_all) ## 390
d=d[assoc_ID %in% studies_in_all]

#5+7+10+41+258+31+16
# 368

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
