test_that("split_by_condition fails", {

  dat <- cbind(matrix(1, nrow = 5, ncol = 10), matrix(2, nrow = 5, ncol = 10))
  colnames(dat) <- paste0("sample", 1:ncol(dat))

  sample_conditions <- c("control", "case", "case", "control", "case")
  names(sample_conditions) <- paste0("test", 1:length(sample_conditions))

  expect_error(split_by_condition(dat, condition_levels = c("control", "case"),
                                  condition_by_sample = sample_conditions))

  names(sample_conditions) <- paste("sample", 1:length(sample_conditions))
  expect_error(split_by_condition(dat[,1:5], condition_levels = c("control", "case")))
  expect_error(split_by_condition(dat[1:5], condition_by_sample = sample_conditions))
  expect_error(split_by_condition(dat, condition_levels = c("control", "case"),
                                  condition_by_sample = sample_conditions))

  sample_conditions <- c(rep("control", 10), rep("case", 10))
  names(sample_conditions) <- paste0("sample", 1:length(sample_conditions))
  expect_equal(split_by_condition(dat, condition_levels = c("control", "case"),
                                  condition_by_sample = sample_conditions),
               list(control = matrix(1, nrow = 5, ncol = 10, dimnames = list(c(), paste0("sample", 1:10))),
                    case = matrix(2, nrow = 5, ncol = 10, dimnames = list(c(), paste0("sample", 11:20)))))
})
