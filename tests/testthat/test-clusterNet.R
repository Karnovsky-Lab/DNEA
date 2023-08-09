test_that("clusterNet", {

  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))
  fake_dat <- matrix(data = 1, nrow = 5, ncol = 5)

  expect_error(clusterNet(dat, tau = 0.3))
  expect_error(clusterNet(dat, tau = 1.2))
  expect_error(clusterNet(dat, max_iterations = -1))
  expect_error(clusterNet(fake_dat))

  expect_s4_class(clusterNet(dat), "DNEAresults")
})
