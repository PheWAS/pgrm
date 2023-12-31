% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgrm.R
\name{get_PGRM}
\alias{get_PGRM}
\title{Get instance of PGRM}
\usage{
get_PGRM(
  ancestry = "all",
  build = "hg19",
  phecode_version = "V1.2",
  unique = TRUE
)
}
\arguments{
\item{ancestry}{A string that indicates the ancestry of the PGRM.
Options EAS, EUR, AFR, SAS, AMR, ALL. Default ALL}

\item{build}{A string indicating the genome reference build. Options hg19, hg38. Default is hg19.}

\item{phecode_version}{A string indicating the phecode version. Currently only
V1.2 is supported, which is the default}

\item{unique}{If TRUE, then rows are uniqued by SNP/phecode (only relevant when ancestry == 'ALL')}
}
\value{
A data.table of the PGRM
}
\description{
This function generates a PGRM copy with specified ancestry, build, and phecode version
}
\details{
This function returns a copy of the PGRM. The function assigns a column 'SNP' as
either SNP_hg19 or SNP_hg38, depending on the build argument. The returned PGRM only includes
associations annotated with the specified ancestry, and the risk allele frequency (RAF) is based
on the specified ancestry.
}
\examples{
library(pgrm)

## Get a copy of the PGRM for build hg19, East Asian ancestry
get_PGRM(build = 'hg19', ancestry = 'EAS')
}
\seealso{
\link{PGRM_ALL}, \code{\link[=annotate_results]{annotate_results()}}
}
