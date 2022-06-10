
get_PGRM(build="hg38x")
get_PGRM(build="hg19")
get_PGRM()
results_MGI

foo=annotate_results(results_MGI,build="hg19",calculate_power=T)
nrow(foo)
head(foo)
?annotate_results

RR=get_RR_power(foo)
RR
