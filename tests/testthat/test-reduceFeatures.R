test_that("reduceFeatures", {

  dat <- readRDS(test_path("testdata", "test-TEDDYresults.rds"))

  TEDDY_groups <- data.frame(features = rownames(expressionData(TEDDYresults, normalized = FALSE)),
                             groups = rownames(expressionData(TEDDYresults, normalized = FALSE)),
                             row.names = rownames(expressionData(TEDDYresults, normalized = FALSE)))

  TEDDY_groups$groups[TEDDY_groups$groups %in% c("isoleucine", "leucine", "valine")] <- "BCAAs"
  TEDDY_groups$groups[grep("acid", TEDDY_groups$groups)] <- "fatty_acids"

  expect_error(reduceFeatures(dat, method = "correlation", correlation_threshold = 0.9))
  expect_error(reduceFeatures(dat, method = "hybrid", correlation_threshold = 0.7))
  expect_error(reduceFeatures(dat, method = "hybrid", correlation_threshold = 0.9, feature_groups = TEDDY_groups))
  expect_error(reduceFeatures(dat, method = "knowledge", correlation_threshold = 0.9))
  expect_warning(reduceFeatures(dat, method = "knowledge", correlation_threshold = 0.7, feature_groups = TEDDY_groups))

  expect_output(reduceFeatures(dat, method = "correlation", correlation_threshold = 0.7))
  expect_output(reduceFeatures(dat, method = "hybrid", correlation_threshold = 0.7, feature_groups = TEDDY_groups))
})
