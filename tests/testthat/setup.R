library('data.table')
library('foreach')
library('qs')
library('BEDMatrix')
registerDoSEQ()

snapshot = function(xObs, path) {
  if (file.exists(path)) {
    xExp = qs::qread(path)
  } else {
    qs::qsave(xObs, path)
    xExp = xObs}
  return(xExp)}

dataDir = 'data'
