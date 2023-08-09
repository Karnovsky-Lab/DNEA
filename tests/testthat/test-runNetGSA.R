test_that("runNetGSA", {

  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))
  fake_dat <- matrix(data = 1, nrow = 5, ncol = 5)

  expect_error(runNetGSA(dat, min_size = -5))
  expect_error(runNetGSA(fake_dat))

  expect_s4_class(runNetGSA(dat), "DNEAresults")

})
