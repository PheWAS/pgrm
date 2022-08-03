#' The phenotype genotype referene map (PGRM)
#'
#' This table provides GWAS associations results from the GWAS catalog, annotated by phecodes, ancestry, etc.
#'
#' @format A data.table with the following columns:
#'
#' * `assoc_ID`: Integer vector of the unique identifier for each association
#' * `SNP_hg19`: Character vector variants encoded in hg19 genome reference build
#' * `SNP_hg38`: Character vector variants encoded in hg38 genome reference build
#' * `rsID`: Character indicating the Reference SNP cluster ID
#' * `risk_allele_dir`: Character vector specifying the risk allele direction
#'   (**alt**: Alternate allele, **ref**: Reference allele)
#' * `risk_allele`: Character specifying the risk allele
#' * `phecode`: Character vector of the phecode for the association
#' * `phecode_string`: Character vector of the phecode string label
#' * `category_string`: Character vector of the phecode category
#' * `ancestry`: Character vector of the ancestery of the source GWAS using 1000
#' Genomes superpopulations ancestry groupings (**AFR**: African, **EAS**: East Asian,
#' **EUR**: European, **AMR**: AdMixed American, **SAS**: South Asian)
#' * `cat_LOG10_P`: Numeric vector of the -log10(P) of the association from the catalog
#' * `cat_OR`: Numeric value of the odds ratio of the association
#' * `cat_L95`: Numeric value of the 95% lower confidence interval of the association
#' * `cat_U95`: Numeric value of the 95% upper confidence interval of the association
#' * `Study_accession`: Character vector of the study accession ID from the GWAS catalog
#' * `pubmedid`: Numeric value of the PubMed ID for source publication of association
#' * `pub_count`: Numeric value indicating how many times association has been published in GWAS catalog
#' * `first_pub_date`: Date of the first date the association was published
#' * `cases_needed_AFR`: Numeric value of the estimated cases needed at 80% power for African ancestry cohort
#' * `cases_needed_EAS`: Numeric value of the estimated cases needed at 80% power for East Asian ancestry cohort
#' * `cases_needed_AMR`: Numeric value  of the estimated cases needed at 80% power for Latino/Admixed American ancestry cohort
#' * `cases_needed_SAS`: Numeric value  of the estimated cases needed at 80% power for South Asian ancestry cohort
#' * `cases_needed_ALL`: Numeric value  of the estimated cases needed at 80% power for multi-ancestry cohort
#' * `AFR_RAF`: Numeric value of the African ancestry risk allele frequency from gnomAD
#' * `EUR_RAF`: Numeric value of the European (non-Finnish) ancestry risk allele frequency from gnomAD
#' * `EAS_RAF`: Numeric value of the East Asian ancestry risk allele frequency from gnomAD
#' * `AMR_RAF`: Numeric value of the Latino/Admixed American ancestry risk allele frequency from gnomAD
#' * `SAS_RAF`: Numeric value of the South Asian ancestry risk allele frequency from gnomAD
#' * `ALL_RAF`: Numeric value of the risk allele frequency from gnomAD
#'
#' @details The odds ratio and 95% confidence intervals (columns cat_OR, cat_L95,
#' and cat_U95) are reported relative to the risk allele which is specified in the
#' column risk_allele.
#'
#' @source Allele frequencies: <https://gnomad.broadinstitute.org/>
#' Summary statistics: <https://www.ebi.ac.uk/gwas/>
#'
#' @seealso [get_PGRM()]
'PGRM_ALL'

#' Information about phecodes
#'
#' The data are artificial and do not correspond to real patients.
#'
#' @format A data table with the following columns:
#'
#' * `phecode`: Character vector of phecode
#' * `phecode_string`: Character vector or string label for phecode
#' * `sex`: Indicates if phecode is sex specific. Values `Female`, `Male`, or `Both`
#' * `leaf`: A boolean value indicating that the phecode is a leaf (has no children)
#' * `phecode_top`: Character vector of the parent phecode
#' * `phecode_top`: Character vector of the category label
#'
#' @source <https://phewascatalog.org>
#'
#' @seealso [get_pheno()]
'phecode_info'

#' Exclude ranges for phecode controls
#'
#' The data are artificial and do not correspond to real patients.
#'
#' @format A data table with the following columns:
#'
#' * `phecode`: Character vector of phecode
#' * `range_start`: Numeric vector of the start of the exclude range
#' * `range_end`: Numeric vector of the end of the exclude range
#'
#' @source <https://phewascatalog.org>
#'
#' @seealso [get_pheno()]
'exclude_ranges'

#' Summary statistics from BioBank Japan (BBJ)
#'
#' This table provides BBJ GWAS summary statistics for SNPs and phecodes that are in the PGRM
#'
#' @format A data.table with the following columns:
#'
#' * `SNP`: String of SNP identifier, format CHR:position:ref:alt. Build = hg19.
#' * `phecode`: String of a phecode
#' * `cases`: Integer of the number of cases for specified phecode
#' * `controls`: Integer of the number of controls for specified phecode
#' * `odds_ratio`: Numeric odds ratio
#' * `P`: Numeric raw P-value
#' * `L95`: Numeric 95% lower confidence interval
#' * `U95`: Numeric 95% upper confidence interval
#' * `cohort_match`: Integer that is 1 if the source of the PGRM association includes BBJ, 0 if not
#'
#' @source <https://pheweb.jp>
#'
'results_BBJ'

#' Summary statistics from UK Biobank (UKBB)
#'
#' This table provides UK Biobank GWAS summary statistics for SNPs and phecodes that are in the PGRM
#'
#' @format A data.table with the following columns:
#'
#' * `SNP`: String of SNP identifier, format CHR:position:ref:alt. Build = hg38
#' * `phecode`: String of a phecode
#' * `cases`: Integer of the number of cases for specified phecode
#' * `controls`: Integer of the number of controls for specified phecode
#' * `odds_ratio`: Numeric odds ratio
#' * `P`: Numeric raw P-value
#' * `L95`: Numeric 95% lower confidence interval
#' * `U95`: Numeric 95% upper confidence interval
#' * `cohort_match`: Integer that is 1 if the source of the PGRM association includes UKBB, 0 if not
#'
#' @source <https://pheweb.org/UKB-TOPMed/>
#'
'results_UKBB'


#' Summary statistics from the Michigan Genomics Initiative (MGI)
#'
#' This table provides MGI GWAS summary statistics for SNPs and phecodes that are in the PGRM
#'
#' @format A data.table with the following columns:
#'
#' * `SNP`: String of SNP identifier, format CHR:position:ref:alt. Build = hg38
#' * `phecode`: String of a phecode
#' * `cases`: Integer of the number of cases for specified phecode
#' * `controls`: Integer of the number of controls for specified phecode
#' * `odds_ratio`: Numeric odds ratio
#' * `P`: Numeric raw P-value
#' * `L95`: Numeric 95% lower confidence interval
#' * `U95`: Numeric 95% upper confidence interval
#' * `cohort_match`: Integer that is 1 if the source of the PGRM association includes MGI, 0 if not
#'
#'
'results_MGI'

#' Summary statistics from the African ancestery cohort in BioVU
#'
#' This table provides BioVU African ancestry GWAS summary statistics for SNPs and phecodes that are in the PGRM
#'
#' @format A data.table with the following columns:
#'
#' * `SNP`: String of SNP identifier, format CHR:position:ref:alt. Build = hg19
#' * `phecode`: String of a phecode
#' * `cases`: Integer of the number of cases for specified phecode
#' * `controls`: Integer of the number of controls for specified phecode
#' * `odds_ratio`: Numeric odds ratio
#' * `P`: Numeric raw P-value
#' * `L95`: Numeric 95% lower confidence interval
#' * `U95`: Numeric 95% upper confidence interval
#' * `cohort_match`: Integer that is 1 if the source of the PGRM association includes BioVU, 0 if not
#'
#'
'results_BioVU_AFR'

#' Summary statistics from the European ancestery cohort in BioVU
#'
#' This table provides BioVU European ancestry GWAS summary statistics for SNPs and phecodes that are in the PGRM
#'
#' @format A data.table with the following columns:
#'
#' * `SNP`: String of SNP identifier, format CHR:position:ref:alt. Build = hg19
#' * `phecode`: String of a phecode
#' * `cases`: Integer of the number of cases for specified phecode
#' * `controls`: Integer of the number of controls for specified phecode
#' * `odds_ratio`: Numeric odds ratio
#' * `P`: Numeric raw P-value
#' * `L95`: Numeric 95% lower confidence interval
#' * `U95`: Numeric 95% upper confidence interval
#' * `cohort_match`: Integer that is 1 if the source of the PGRM association includes BioVU, 0 if not
#'
#'
'results_BioVU_EUR'

#' Summary statistics from all five test cohorts with annotations
#'
#' This table summary statistics for PGRM association in BioVU (EUR and AFR), BBJ, MGI, and UKBB
#'
#' @format A data.table with the following columns:
#'
#' * `assoc_ID`: Integer vector of the unique identifier for each association
#' * `phecode`: Character vector of the phecode for the association
#' * `phecode_string`: Character vector of the phecode string label
#' * `category_string`: Character vector of the phecode category
#' * `cases`: Integer of the number of cases for specified phecode/cohort
#' * `controls`: Integer of the number of controls for specified phecode/cohort
#' * `odds_ratio`: Numeric odds ratio for test cohort association
#' * `P`: Numeric raw P-value for test cohort association
#' * `L95`: Numeric 95% lower confidence interval for test cohort association
#' * `U95`: Numeric 95% upper confidence interval for test cohort association
#' * `ancestry`: Character vector of the ancestery of the source GWAS using 1000
#' Genomes superpopulations ancestry groupings (**AFR**: African, **EAS**: East Asian,
#' **EUR**: European, **AMR**: AdMixed American, **SAS**: South Asian)
#' * `rsID`: Character indicating the Reference SNP cluster ID
#' * `risk_allele_dir`: Character indicating the specifying the risk allele direction
#'   (**alt**: Alternate allele, **ref**: Reference allele)
#' * `RAF`: GnomAD frequency of risk allele, matched to ancestry
#' * `cat_LOG10_P`: Numeric vector of the -log10(P) of the association from the catalog
#' * `cat_OR`: Numeric vector of the odds ratio of the association
#' * `cat_L95`: Numeric vector of the 95% lower confidence interval of the association
#' * `cat_U95`: Numeric vector of the 95% upper confidence interval of the association
#' * `Study_accession`: Character vector of the study accession ID from the GWAS catalog
#' * `Powered`: Boolean value indicating if the association is powered at >=80%
#' * `rep`: Boolean value indicating if the association replicated
#' * `CI_overlap`: Character vector indicating the overlap of the confidence intervals from
#'   the catalog and test association
#' * `CI_overlap`: Character vector indicating the test cohort
#' * `pub_count`: The number of times the association has been reported in the catalog
#' * `first_pub_date`: The first date the association was published
#'
'benchmark_results'
