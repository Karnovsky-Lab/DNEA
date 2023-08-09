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
