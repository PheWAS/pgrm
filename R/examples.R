example1 = function() {
  ex = "
  @examples
  # get an instance of PGRM with build hg19 SNPs and East Asian ancestry
  PGRM_EAS = get_PGRM(build='hg19',ancestry='EAS')
  "
  return(strsplit(ex, split = '\n')[[1L]])}
