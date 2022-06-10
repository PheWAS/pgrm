checkBuild = function(build) {
  assert(build %in% c('hg19','hg37'))
  #assertDataTable(demos)
  #assert(anyDuplicated(demos$person_id) == 0)
  invisible()
}
