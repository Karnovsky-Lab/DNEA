test_that("includeMetadata", {

  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))
  metadat <- readRDS(test_path("testdata", "test-T1Dmeta.rds"))

  metadat2 <- matrix(data = 1, nrow = 5, ncol = 10, dimnames = list(paste0("row", seq(1,5)), paste0("col", seq(1,10))))

  expect_error(includeMetadata(dat, type = "features", metadata = metadat))
  expect_error(includeMetadata(dat, type = "samples", metadata = metadat2))

  expect_s4_class(includeMetadata(dat, type = "samples", metadata = metadat), "DNEA")
})

test_that("getNetworkFiles", {

  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))
  expect_error(withr::with_tempdir(getNetworkFiles(dat, file_path = 7), clean = TRUE))
  expect_no_condition(withr::with_tempdir(getNetworkFiles(dat), clean = TRUE))


})

test_that("plotNetworks", {

  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))

  expect_error(plotNetworks(dat, type = "test"))
  expect_error(plotNetworks(dat, type = "group_networks", subtype = 1))
  expect_error(plotNetworks(dat, type = "sub_networks", subtype = "All"))

  expect_no_condition(withr::with_tempdir(plotNetworks(dat, type = "group_networks", subtype = "All"), clean = TRUE))
  expect_no_condition(withr::with_tempdir(plotNetworks(dat, type = "sub_networks", subtype = 1), clean = TRUE))
})

test_that("filterNetworks", {

  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))

  expect_error(filterNetworks(dat, pcor = 0.5, top_percent_edges = 0.5))
  expect_error(filterNetworks(dat, pcor = 1.2))
  expect_error(filterNetworks(dat, pcor = -0.3))
  expect_error(filterNetworks(dat, top_percent_edges = 1.2))
  expect_error(filterNetworks(dat, top_percent_edges = -0.3))
  expect_s4_class(filterNetworks(dat, pcor = 0.3), "DNEA")
})
