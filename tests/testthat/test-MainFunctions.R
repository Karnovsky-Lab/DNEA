test_that("BICtune", {

  library(BiocParallel)
  BP_plan <- BiocParallel::SerialParam(RNGseed = 417)
  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))

  expect_error(BICtune(dat, lambda_values = list(0, 0.5, 0.2, 1.2), BPPARAM = BP_plan))
  expect_error(BICtune(dat, lambda_values = list(0, 0.5, 0.2, -0.2), BPPARAM = BP_plan))

  expect_output(BICtune(dat, lambda_values = list(0.0391997757765919), BPPARAM = BP_plan))
})

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

test_that("getNetworks", {

  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))
  fake_dat <- matrix(data = 1, nrow = 5, ncol = 5)

  expect_error(getNetworks(dat, optimal_lambda = -0.5))
  expect_error(getNetworks(dat, optimal_lambda = 1.2))
  expect_error(getNetworks(dat, optimal_lambda = -0.5))
  expect_error(getNetworks(dat, eps_threshold = -0.5))
  expect_error(getNetworks(dat, eps_threshold = 1.2))
  expect_error(getNetworks(fake_dat))

  expect_output(getNetworks(dat))
})

test_that("clusterNet", {

  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))
  fake_dat <- matrix(data = 1, nrow = 5, ncol = 5)

  expect_error(clusterNet(dat, tau = 0.3))
  expect_error(clusterNet(dat, tau = 1.2))
  expect_error(clusterNet(dat, max_iterations = -1))
  expect_error(clusterNet(fake_dat))

  expect_output(clusterNet(dat))
})
