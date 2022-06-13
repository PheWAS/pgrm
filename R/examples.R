example_get_PGRM = function() {
  ex = "
  @examples
  library('pgrm')

  # map ICD codes to phecodes
  phecodeOccurrences = getPhecodeOccurrences(icdSample)

  # calculate weights
  weights = getWeights(demoSample, phecodeOccurrences)

  # OMIM disease IDs for which to calculate phenotype risk scores
  diseaseId = 154700

  # map diseases to phecodes
  diseasePhecodeMap = mapDiseaseToPhecode()

  # calculate scores
  scores = getScores(
  demoSample, phecodeOccurrences, weights, diseasePhecodeMap[disease_id == diseaseId])

  # calculate residual scores
  rscores = getResidualScores(demoSample, scores, lmFormula = ~ sex)
  "
  return(strsplit(ex, split = '\n')[[1L]])}
