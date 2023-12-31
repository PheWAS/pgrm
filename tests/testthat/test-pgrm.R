test_that('get_PGRM with build hg38', {
  PGRM_observed = get_PGRM(build = 'hg38')
  PGRM_expected = snapshot(PGRM_observed, file.path('snapshots', 'PGRM_hg38.qs'))
  expect_equal(PGRM_observed, PGRM_expected)
})


test_that('test annotate_results with MGI data', {
  MGI_annotate_observed = annotate_results(
    head(results_MGI, n = 20), build = 'hg38', ancestry = 'EUR', calculate_power = TRUE, LOUD = FALSE)
  MGI_annotate_expected = snapshot(MGI_annotate_observed, file.path('snapshots', 'MGI_annotate.qs'))
  expect_equal(MGI_annotate_observed, MGI_annotate_expected)

  MGI_RR_observed = get_RR(MGI_annotate_observed, LOUD = FALSE, include = 'powered')
  MGI_RR_expected = snapshot(MGI_RR_observed, file.path('snapshots', 'MGI_RR_powered.qs'))
  expect_equal(MGI_RR_observed, MGI_RR_expected)

  MGI_RR_ALL_observed = get_RR(MGI_annotate_observed, LOUD = FALSE, include = 'all')
  MGI_RR_ALL_expected = snapshot(MGI_RR_ALL_observed, file.path('snapshots', 'MGI_RR_all.qs'))
  expect_equal(MGI_RR_ALL_observed, MGI_RR_ALL_expected)

  MGI_AE_observed = get_AER(MGI_annotate_observed, LOUD = FALSE)
  MGI_AE_expected = snapshot(MGI_AE_observed, file.path('snapshots', 'MGI_AE.qs'))
  expect_equal(MGI_AE_observed, MGI_AE_expected)

  MGI_powered_observed = get_powered_rate(MGI_annotate_observed)
  MGI_powered_expected = snapshot(MGI_AE_observed, file.path('snapshots', 'MGI_powered.qs'))
})

