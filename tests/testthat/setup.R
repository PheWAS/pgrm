library('data.table')
#registerDoSEQ()

snapshot = function(xObs, path) {
  if (file.exists(path)) {
    xExp = qs::qread(path)
  } else {
    qs::qsave(xObs, path)
    xExp = xObs}
  return(xExp)}

dataDir = 'data'

phecode_table_test = data.table(
  person_id = c(1, rep(2L, 3), 3, rep(4L, 3)),
  phecode = c('250.2', '250.2', '174.11', '153', '250', '250.1', '427.21', '185'),
  N= c(2,2,8,7,8,5,3,7)
)

demos_table_test = data.table(
  person_id = 1:4, sex = c('F', 'M', 'F', 'M'),
  last_age = c(28772,18028,11636,14589))
