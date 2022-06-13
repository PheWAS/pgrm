
library(pgrm)
get_PGRM(build="hg38")

get_PGRM(build="hg19x")
get_PGRM()


foo=annotate_results(head(results_MGI,n=20),build="hg38",calculate_power=T,LOUD=F)
nrow(foo)
head(foo)

?annotate_results

RR=get_RR_power(foo)
RR
