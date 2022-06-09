#' The phenotype genotype referene map (PGRM)
#'
#' This table provides GWAS associations results from the GWAS catalog, annotated by phecodes, ancestry, etc.
#'
#' @format A data.table with the following columns:
#'
#' * `assoc_ID`: Integer vector of the unique identifier for each association
#' * `SNP_hg19`: Character vector variants encoded in hg19 genome reference build
#' * `SNP_hg38`: Character vector variants encoded in hg38 genome reference build
#' * `ancestry`: Character vector of the ancestery of the source GWAS using 1000 Genomes superpopulations ancestry groupings
#'   (**AFR**: African, **EAS**: East Asian, **EUR**: European, **AMR**: AdMixed American, **SAS**: South Asian)
#' * `rsID`: Character indicating the Reference SNP cluster ID
#' * `risk_allele_dir`: Character indicating the specifying the risk allele
#'   (**alt**: Alternate allele, **ref**: Reference allele)
#' * `phecode`: Character vector of the phecode for the association
#' * `phecode_string`: Character vector of the phecode string label
#' * `best_log_P`: Numeric vector of the -log10(P) of the association
#' * `cat_OR`: Numeric vector of the odds ratio of the association
#' * `cat_L95`: Numeric vector of the 95% lower confidence interval of the association
#' * `cat_U95`: Numeric vector of the 95% upper confidence interval of the association
#' * `Study_accession`: Character vector of the study accession ID from the GWAS catalog
#' * `cases_needed_AFR`: Integer vector of the estimated cases needed at 80% power for African ancestry cohort
#' * `cases_needed_EAS`: Integer vector of the estimated cases needed at 80% power for East Asian ancestry cohort
#' * `cases_needed_AMR`: Integer vector of the estimated cases needed at 80% power for AdMixed American ancestry cohort
#' * `cases_needed_SAS`: Integer vector of the estimated cases needed at 80% power for South Asian ancestry cohort
#' * `cases_needed_ALL`: Integer vector of the estimated cases needed at 80% power for multi-ancestry cohort
#' * `category_string`: Character vector of the phecode category
#'
#' @seealso [get_PGRM()]]
'PGRM_ALL'

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
#' @seealso [get_pheno()]]
'exclude_ranges'

