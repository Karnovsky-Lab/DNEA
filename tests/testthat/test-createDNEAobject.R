test_that("createDNEAobject", {

  TEDDY <- readRDS(test_path("testdata", "test-TEDDY.rds"))
  T1Dmeta <- readRDS(test_path("testdata", "test-T1Dmeta.rds"))

  unnamed_labels <- factor(T1Dmeta$group, levels = c("DM:control", "DM:case"))
  group_labels <- factor(T1Dmeta$group, levels = c("DM:control", "DM:case"))
  names(group_labels) <- rownames(T1Dmeta)

  character_dat <- matrix(data = "1", nrow = nrow(TEDDY), ncol = ncol(TEDDY), dimnames = list(rownames(TEDDY), colnames(TEDDY)))
  unnamed_dat <- matrix(data = "1", nrow = nrow(TEDDY), ncol = ncol(TEDDY))

  expect_error(createDNEAobject(project_name = 5, expression_data = TEDDY, group_labels = group_labels))
  expect_error(createDNEAobject(project_name = "testing", expression_data = TEDDY, group_labels = unnamed_labels))
  expect_error(createDNEAobject(project_name = "testing", expression_data = character_dat, group_labels = group_labels))
  expect_error(createDNEAobject(project_name = "testing", expression_data = unnamed_dat, group_labels = group_labels))

})
