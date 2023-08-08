test_that("includeMetadata", {

  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))
  metadat <- readRDS(test_path("testdata", "test-T1Dmeta.rds"))

  metadat2 <- matrix(data = 1, nrow = 5, ncol = 10, dimnames = list(paste0("row", seq(1,5)), paste0("col", seq(1,10))))

  expect_error(includeMetadata(dat, type = "feature", metadata = metadat))
  expect_error(includeMetadata(dat, type = "sample", metadata = metadat2))

  expect_s4_class(includeMetadata(dat, type = "sample", metadata = metadat), "DNEAresults")
})

test_that("getNetworkFiles", {

  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))
  expect_no_condition(getNetworkFiles(dat, file_path = test_path()))

  withr::defer(unlink(paste0(test_path(),"/.", dat@project_name,'_edgelist_',Sys.Date(),'.csv')))
  withr::defer(unlink(paste0(test_path(),"/.", dat@project_name,'_edgelist_',Sys.Date(),'.csv')))
})

test_that("plotNetworks", {

  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))

  expect_error(plotNetworks(dat, type = "test"))
  expect_error(plotNetworks(dat, type = "group_networks", subtype = 1))
  expect_error(plotNetworks(dat, type = "subnetworks", subtype = "All"))

  expect_no_condition(plotNetworks(dat, type = "group_networks", subtype = "All"))
  expect_no_condition(plotNetworks(dat, type = "subnetworks", subtype = 1))

  withr::defer(unlink(paste0(test_path(),"/Rplots.pdf")))
})

test_that("filterNetworks", {

  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))

  expect_error(filterNetworks(dat, pcor = 0.5, top_percent_edges = 0.5))
  expect_error(filterNetworks(dat, pcor = 1.2))
  expect_error(filterNetworks(dat, pcor = -0.3))
  expect_error(filterNetworks(dat, top_percent_edges = 1.2))
  expect_error(filterNetworks(dat, top_percent_edges = -0.3))
  expect_s4_class(filterNetworks(dat, pcor = 0.3), "DNEAresults")
})





