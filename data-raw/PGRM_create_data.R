library('data.table')

rawDir = 'data-raw'

########################
# PGRM. Open the map and allele frequencies. Merge on SNP_hg19

 PGRM = fread(
   file.path(rawDir, 'PGRM.csv'),
   colClasses = list(character = 'phecode'),
   na.strings="-")
setkeyv(PGRM, c('assoc_ID', 'phecode','SNP_hg19'))
 PGRM_AF = fread(file.path(rawDir, 'PGRM_AF.csv'))
  setkeyv(PGRM_AF, c('SNP_hg19'))
 PGRM=merge(PGRM,PGRM_AF,by="SNP_hg19")

 usethis::use_data(PGRM, overwrite = TRUE)

#######################
# phecode exclude ranges for PheWAS

 exclude_ranges = fread(
   file.path(rawDir, 'phecode_ranges.csv'),
   colClasses = list(character = 'phecode'))
 setkeyv(exclude_ranges, c( 'phecode'))

 usethis::use_data(exclude_ranges, overwrite = TRUE)

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
# Results file - MGI

 results_UKBB = fread(
   file.path(rawDir, 'results_UKBB.csv'),
   colClasses = list(character = 'phecode'))
 setkeyv(results_UKBB, c('SNP', 'phecode'))

 usethis::use_data(results_UKBB, overwrite = TRUE)
