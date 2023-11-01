test_that("BICtune", {

  library(BiocParallel)
  BP_plan <- BiocParallel::SerialParam(RNGseed = 417)
  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))
  fakedat <- matrix(1, nrow = 10, ncol = 10)

  expect_error(clusterNet(object = object, tau = 0.3))
  expect_error(clusterNet(object = object, tau = 1.2))
  expect_error(clusterNet(object = object, tau = 0.5, max_iterations = -5))
  expect_error(clusterNet(object = object, tau = 0.5, max_iterations = 5, verbose = 7))
  expect_error(clusterNet(object = fakedat))
})
