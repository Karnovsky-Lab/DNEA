test_that("split_by_condition fails", {

  dat <- cbind(matrix(1, nrow = 5, ncol = 10), matrix(2, nrow = 5, ncol = 10))
  expect_error(split_by_condition(dat, condition_levels = c("control", "case")))
  expect_error(split_by_condition(dat, condition_by_sample = c("control", "case", "case", "control", "case")))
  expect_error(split_by_condition(dat, condition_levels = c("control", "case"),
                                  condition_by_sample = c("control", "case", "case", "control", "case")))

  expect_equal(split_by_condition(dat, condition_levels = c("control", "case"),
                                  condition_by_sample = c(rep("control", 10), rep("case", 10))),
               list(control = matrix(1, nrow = 5, ncol = 10), case = matrix(2, nrow = 5, ncol = 10)))
})

