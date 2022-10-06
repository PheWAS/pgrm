library('data.table')

rawDir = 'data-raw'

########################
# The PGRM. Add allele frequnecies

pgrm_all = fread(
  file.path(rawDir, 'PGRM_ALL.csv'),
  colClasses = list(character = 'phecode'),
  na.strings = '-')
setkeyv(pgrm_all, c('assoc_ID', 'phecode', 'SNP_hg19','SNP_hg38'))

PGRM_AF = fread(file.path(rawDir, 'PGRM_AF.csv'))
setkeyv(PGRM_AF, c('SNP_hg19'))

pgrm_all = merge(pgrm_all, PGRM_AF, by = 'SNP_hg19')
pgrm_all[phecode == '555', phecode_string := 'Inflammatory bowel disease']

PGRM_pubinfo = fread(
  file.path(rawDir, 'PGRM_pubinfo.csv'))
setkeyv(PGRM_pubinfo, c('assoc_ID'))
pgrm_all = merge(pgrm_all, PGRM_pubinfo, by = 'assoc_ID')

pgrm_all$pub_date=as.Date(pgrm_all$pub_date)
nrow(pgrm_all)
head(pgrm_all)
ncol(pgrm_all)
pgrm_all[risk_allele_dir=="ref"]$AFR_RAF=1-pgrm_all[risk_allele_dir=="ref"]$AFR_RAF
pgrm_all[risk_allele_dir=="ref"]$EUR_RAF=1-pgrm_all[risk_allele_dir=="ref"]$EUR_RAF
pgrm_all[risk_allele_dir=="ref"]$EAS_RAF=1-pgrm_all[risk_allele_dir=="ref"]$EAS_RAF
pgrm_all[risk_allele_dir=="ref"]$AMR_RAF=1-pgrm_all[risk_allele_dir=="ref"]$AMR_RAF
pgrm_all[risk_allele_dir=="ref"]$SAS_RAF=1-pgrm_all[risk_allele_dir=="ref"]$SAS_RAF
pgrm_all[risk_allele_dir=="ref"]$ALL_RAF=1-pgrm_all[risk_allele_dir=="ref"]$ALL_RAF
head(pgrm_all)
pgrm_all$first_pub_date=NULL
pgrm_all[, c('cases_needed_AFR','cases_needed_EAS','cases_needed_EUR','cases_needed_AMR','cases_needed_SAS','cases_needed_ALL'):=NULL]
setcolorder(pgrm_all, c('assoc_ID', 'SNP_hg19', 'SNP_hg38','rsID', 'risk_allele_dir','risk_allele',
                  'phecode','phecode_string','category_string','ancestry',
                  'cat_LOG10_P','cat_OR','cat_L95','cat_U95','Study_accession','pubmedid','pub_count','pub_date',
                  'AFR_RAF','EUR_RAF','EAS_RAF','AMR_RAF','SAS_RAF','ALL_RAF'))

setnames(pgrm_all, 'assoc_ID', 'assoc_id')
setnames(pgrm_all, 'SNP_hg19', 'snp_hg19')
setnames(pgrm_all, 'SNP_hg38', 'snp_hg38')
setnames(pgrm_all, 'rsID', 'rsid')
setnames(pgrm_all, 'cat_LOG10_P', 'cat_log10_p')
setnames(pgrm_all, 'cat_OR', 'cat_or')
setnames(pgrm_all, 'cat_L95', 'cat_l95')
setnames(pgrm_all, 'cat_U95', 'cat_u95')
setnames(pgrm_all, 'Study_accession', 'study_accession')
setnames(pgrm_all, 'AFR_RAF', 'afr_raf')
setnames(pgrm_all, 'EUR_RAF', 'eur_raf')
setnames(pgrm_all, 'EAS_RAF', 'eas_raf')
setnames(pgrm_all, 'AMR_RAF', 'amr_raf')
setnames(pgrm_all, 'SAS_RAF', 'sas_raf')
setnames(pgrm_all, 'ALL_RAF', 'all_raf')
pgrm_all$ancestry = tolower(pgrm_all$ancestry)
usethis::use_data(pgrm_all, overwrite = TRUE)

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
# Raw results files

results_BioVU_EUR = fread(
     file.path(rawDir, 'results_BioVU_EUR.csv'),
     colClasses = list(character = 'phecode'))
 results_BioVU_EUR$dataset = 'biovu_eur'
 setnames(results_BioVU_EUR,'SNP','SNP_hg19')
 results_BioVU_EUR = merge(results_BioVU_EUR, unique(PGRM_ALL[,c('SNP_hg19','SNP_hg38')]),by='SNP_hg19')
 setnames(results_BioVU_EUR,'SNP_hg38','SNP')
 results_BioVU_EUR$SNP_hg19 = NULL
 setcolorder(results_BioVU_EUR, c('SNP','phecode','cases','controls','odds_ratio','P','L95',
                         'U95','dataset','cohort_match'))

results_BioVU_AFR = fread(
    file.path(rawDir, 'results_BioVU_AFR.csv'),
    colClasses = list(character = 'phecode'))
 results_BioVU_AFR$dataset = 'biovu_afr'
 setnames(results_BioVU_AFR,'SNP','SNP_hg19')
 results_BioVU_AFR = merge(results_BioVU_AFR, unique(PGRM_ALL[,c('SNP_hg19','SNP_hg38')]),by='SNP_hg19')
 setnames(results_BioVU_AFR,'SNP_hg38','SNP')
 results_BioVU_AFR$SNP_hg19 = NULL
 setcolorder(results_BioVU_AFR, c('SNP','phecode','cases','controls','odds_ratio','P','L95',
                                  'U95','dataset','cohort_match'))

results_BBJ = fread(
    file.path(rawDir, 'results_BBJ.csv'),
    colClasses = list(character = 'phecode'))
 results_BBJ$dataset = 'bbj'
 setnames(results_BBJ,'SNP','SNP_hg19')
 results_BBJ = merge(results_BBJ, unique(PGRM_ALL[,c('SNP_hg19','SNP_hg38')]),by='SNP_hg19')
 setnames(results_BBJ,'SNP_hg38','SNP')
 results_BBJ$SNP_hg19 = NULL
 setcolorder(results_BBJ, c('SNP','phecode','cases','controls','odds_ratio','P','L95',
                                 'U95','dataset','cohort_match'))

results_MGI = fread(
    file.path(rawDir, 'results_MGI.csv'),
    colClasses = list(character = 'phecode'))
results_MGI$dataset = 'mgi'
setcolorder(results_MGI, c('SNP','phecode','cases','controls','odds_ratio','P','L95',
                           'U95','dataset','cohort_match'))

results_UKBB = fread(
    file.path(rawDir, 'results_UKBB.csv'),
    colClasses = list(character = 'phecode'))
results_UKBB$dataset = 'ukb'
setcolorder(results_UKBB, c('SNP','phecode','cases','controls','odds_ratio','P','L95',
                           'U95','dataset','cohort_match'))

summary_stats = rbind(results_BioVU_EUR, results_BioVU_AFR)
summary_stats = rbind(summary_stats, results_BBJ)
summary_stats = rbind(summary_stats, results_MGI)
summary_stats = rbind(summary_stats, results_UKBB)
setnames(summary_stats, 'SNP','snp')
setnames(summary_stats, 'P','p')
setnames(summary_stats, 'L95','l95')
setnames(summary_stats, 'U95','u95')

setkeyv(summary_stats, c('dataset','snp','phecode'))
usethis::use_data(summary_stats, overwrite = TRUE)

#######################
# Benchmark results

summary_stats_anno = fread(
   file.path(rawDir, 'annotated_results.csv'),
   colClasses = list(character = 'phecode'))
 setkeyv(summary_stats_anno, c('assoc_ID'))


 summary_stats_anno$first_pub_date=as.Date(summary_stats_anno$first_pub_date)
 usethis::use_data(summary_stats_anno, overwrite = TRUE)

 #######################
 # ICD to phecode map (version 1.2)

 icdPhecodeMap_V1_2 = fread(
   file.path(rawDir, 'icd_to_phecode_V1_2.csv'),
   colClasses = list(character = 'phecode','icd'))
 setkeyv(icdPhecodeMap_V1_2, c('icd','flag','phecode'))

 usethis::use_data(icdPhecodeMap_V1_2, overwrite = TRUE)


 #######################
 # Example ICD file

 icdExampleTable = fread(
   file.path(rawDir, 'icds_example.csv'),
   colClasses = list(character = 'icd'))
 setkeyv(icdExampleTable, c('icd','flag'))
 icdExampleTable$V1 = NULL
 usethis::use_data(icdExampleTable, overwrite = TRUE)
