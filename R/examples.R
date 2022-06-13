example1 = function() {
  ex = "
  @examples
  # get an instance of PGRM with build hg19 SNPs and East Asian ancestry
  PGRM_EAS = get_PGRM(build='hg19',ancestry='EAS')
  "
  return(strsplit(ex, split = '\n')[[1L]])
}

example2 = function() {
  return("foobar")
}

example3 = function() {
  return("foobar2")
}

example4 = function() {
  return("foobar4")
}
