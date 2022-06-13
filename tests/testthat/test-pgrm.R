test_that('get_PGRM with build hg38', {
  # test stuff
  # expect = get_PGRM() ##
  PGRM_observed = get_PGRM(build="hg38")
  PGRM_expected = snapshot(PGRM_observed,"snapshots/PGRM_hg38.qs")
  expect_equal(PGRM_observed,PGRM_expected)

})


test_that('test annotate_results with MGI data', {
  MGI_annotate_observed=annotate_results(head(results_MGI,n=20),build="hg38",calculate_power=T,LOUD=F)
  MGI_annotate_expected=snapshot(MGI_annotate_observed,"snapshots/MGI_annotate.qs" )
  expect_equal(MGI_annotate_observed,MGI_annotate_expected)
})

