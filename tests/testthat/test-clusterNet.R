test_that("BICtune", {

  library(BiocParallel)
  BP_plan <- BiocParallel::SerialParam(RNGseed = 417)
  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))

  expect_error(BICtune(dat, lambda_values = list(0, 0.5, 0.2, 1.2), BPPARAM = BP_plan))
  expect_error(BICtune(dat, lambda_values = list(0, 0.5, 0.2, -0.2), BPPARAM = BP_plan))

  expect_output(BICtune(dat, lambda_values = list(0.0391997757765919), BPPARAM = BP_plan))
})
