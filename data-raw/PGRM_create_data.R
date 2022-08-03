library('data.table')

rawDir = 'data-raw'

########################
# The PGRM. Add allele frequnecies

PGRM_ALL = fread(
  file.path(rawDir, 'PGRM_ALL.csv'),
  colClasses = list(character = 'phecode'),
  na.strings = '-')
setkeyv(PGRM_ALL, c('assoc_ID', 'phecode', 'SNP_hg19','SNP_hg38'))

PGRM_AF = fread(file.path(rawDir, 'PGRM_AF.csv'))
setkeyv(PGRM_AF, c('SNP_hg19'))

PGRM_ALL = merge(PGRM_ALL, PGRM_AF, by = 'SNP_hg19')
PGRM_ALL[phecode == '555', phecode_string := 'Inflammatory bowel disease']

PGRM_pubinfo = fread(
  file.path(rawDir, 'PGRM_pubinfo.csv'))
setkeyv(PGRM_pubinfo, c('assoc_ID'))
PGRM_ALL = merge(PGRM_ALL, PGRM_pubinfo, by = 'assoc_ID')
PGRM_ALL$first_pub_date=as.Date(PGRM_ALL$first_pub_date)
nrow(PGRM_ALL)
head(PGRM_ALL)
ncol(PGRM_ALL)
PGRM_ALL[risk_allele_dir=="ref"]$AFR_RAF=1-PGRM_ALL[risk_allele_dir=="ref"]$AFR_RAF
PGRM_ALL[risk_allele_dir=="ref"]$EUR_RAF=1-PGRM_ALL[risk_allele_dir=="ref"]$EUR_RAF
PGRM_ALL[risk_allele_dir=="ref"]$EAS_RAF=1-PGRM_ALL[risk_allele_dir=="ref"]$EAS_RAF
PGRM_ALL[risk_allele_dir=="ref"]$AMR_RAF=1-PGRM_ALL[risk_allele_dir=="ref"]$AMR_RAF
PGRM_ALL[risk_allele_dir=="ref"]$SAS_RAF=1-PGRM_ALL[risk_allele_dir=="ref"]$SAS_RAF
PGRM_ALL[risk_allele_dir=="ref"]$ALL_RAF=1-PGRM_ALL[risk_allele_dir=="ref"]$ALL_RAF
head(PGRM_ALL)

setcolorder(PGRM_ALL, c('assoc_ID', 'SNP_hg19', 'SNP_hg38','rsID', 'risk_allele_dir','risk_allele',
                  'phecode','phecode_string','category_string','ancestry',
                  'cat_LOG10_P','cat_OR','cat_L95','cat_U95','Study_accession','pubmedid','pub_count','first_pub_date',
                  'cases_needed_AFR','cases_needed_EAS','cases_needed_EUR','cases_needed_AMR','cases_needed_SAS','cases_needed_ALL',
                  'AFR_RAF','EUR_RAF','EAS_RAF','AMR_RAF','SAS_RAF','ALL_RAF'))

head(PGRM_ALL)

usethis::use_data(PGRM_ALL, overwrite = TRUE)

########################
# PGRM allele frequencies

# PGRM_AF = fread(
#   file.path(rawDir, 'PGRM_AF.csv'))

 #setkeyv(PGRM_AF, c('SNP_hg19'))
 #usethis::use_data(PGRM_AF, overwrite = TRUE)



#######################
# phecode exclude ranges for PheWAS

 exclude_ranges = fread(
   file.path(rawDir, 'phecode_ranges.csv'),
   colClasses = list(character = 'phecode'))
 setkeyv(exclude_ranges, c('phecode'))

 usethis::use_data(exclude_ranges, overwrite = TRUE)

 #######################
 # phecode information including phecode_string, category, whether sex specific

 phecode_info = fread(
   file.path(rawDir, 'phecode_info.csv'),
   colClasses = list(character = 'phecode', character="phecode_top"))
 setkeyv(phecode_info, c( 'phecode'))

 usethis::use_data(phecode_info, overwrite = TRUE)


#######################
# Results file - BioVU EUR

 results_BioVU_EUR = fread(
   file.path(rawDir, 'results_BioVU_EUR.csv'),
   colClasses = list(character = 'phecode'))
 setkeyv(results_BioVU_EUR, c('SNP', 'phecode'))

 usethis::use_data(results_BioVU_EUR, overwrite = TRUE)

#######################
# Results file - BioVU AFR

 results_BioVU_AFR = fread(
   file.path(rawDir, 'results_BioVU_AFR.csv'),
   colClasses = list(character = 'phecode'))
 setkeyv(results_BioVU_AFR, c('SNP', 'phecode'))

 usethis::use_data(results_BioVU_AFR, overwrite = TRUE)

#######################
# Results file - BBJ

 results_BBJ = fread(
   file.path(rawDir, 'results_BBJ.csv'),
   colClasses = list(character = 'phecode'))
 setkeyv(results_BBJ, c('SNP', 'phecode'))

 usethis::use_data(results_BBJ, overwrite = TRUE)

#######################
# Results file - MGI

 results_MGI = fread(
   file.path(rawDir, 'results_MGI.csv'),
   colClasses = list(character = 'phecode'))
 setkeyv(results_MGI, c('SNP', 'phecode'))

 usethis::use_data(results_MGI, overwrite = TRUE)

#######################
# Results file - UKBB

 results_UKBB = fread(
   file.path(rawDir, 'results_UKBB.csv'),
   colClasses = list(character = 'phecode'))
 setkeyv(results_UKBB, c('SNP', 'phecode'))

 usethis::use_data(results_UKBB, overwrite = TRUE)

#######################
# Benchmark results

 benchmark_results = fread(
   file.path(rawDir, 'annotated_results.csv'),
   colClasses = list(character = 'phecode'))
 setkeyv(benchmark_results, c('assoc_ID'))
 head(benchmark_results)

 benchmark_results$first_pub_date=as.Date(benchmark_results$first_pub_date)
 usethis::use_data(benchmark_results, overwrite = TRUE)
