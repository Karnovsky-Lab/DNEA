test_that("stabilitySelection", {

  library(BiocParallel)
  BP_plan <- BiocParallel::SerialParam(RNGseed = 417)
  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))
  fake_dat <- matrix(data = 1, nrow = 5, ncol = 5)

  expect_error(stabilitySelection(dat, subSample = FALSE, nreps = -5, BPPARAM = BP_plan))
  expect_error(stabilitySelection(dat, subSample = "TEST", nreps = -5, BPPARAM = BP_plan))
  expect_error(stabilitySelection(fake_dat, subSample = FALSE, nreps = 5, BPPARAM = BP_plan))
  expect_error(stabilitySelection(dat, subSample = FALSE, optimized_lambda = -0.5, nreps = 5, BPPARAM = BP_plan))

  expect_output(stabilitySelection(dat, subSample = TRUE, nreps = 1, BPPARAM = BP_plan))
})
