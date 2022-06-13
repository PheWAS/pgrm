checkBuild = function(build) {
  if (!build %in% c('hg19','hg38')) {
    stop('Build must equal hg19 or hg38.')
  }
  invisible()}

checkAncestry = function(ancestry) {
  if (!ancestry %in% c('EAS', 'EUR', 'AFR', 'SAS', 'AMR', 'ALL')) {
    stop('Ancestry must be EAS, EUR, AFR, SAS, AMR, or ALL.')
  }
  invisible()}
