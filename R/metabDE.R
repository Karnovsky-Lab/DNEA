#' Perform differential expression analysis using students T-test
#'
#' Performs differential expression on features of the input data provided using student's T-Test and multiple hypothesis
#' correction via Benjamini-Hochberg
#'
#' @param mat An m x n expression matrix with features as columns
#' @param condition_values a vector of the two possible condition values
#' @param conditions a vector containing the condition labels corresponding to the samples in mat
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{createDNEAobject}}
#'
#' @keywords internal
#' @noRd
metabDE <- function(mat,
                    condition_values,
                    conditions){

  ##set necessary parameters and output
  num_features <- ncol(mat)
  num_samples <- nrow(mat)
  feature_info <- data.frame('Features' = colnames(mat))

  ##split data by condition
  cond_data <- split_by_condition(dat = mat,
                                  condition_levels = condition_values,
                                  condition_by_sample = conditions)

  ##perform Differential expression
  assign(paste0('feature_info$foldchange-',condition_values[2], '/', condition_values[1]), rowMeans(cond_data[[2]]) - rowMeans(cond_data[[1]]))
  feature_info$fcdirection <- sapply(1:num_features, function(i) ifelse(get(paste0('feature_info$foldchange-',condition_values[2], '/', condition_values[1]))[i] > 0, "Up", "Down"))
  feature_info$t_statistic <- sapply(1:num_features, function(i) t.test(cond_data[[2]][i, ], cond_data[[1]][i, ], var.equal = FALSE)$statistic)
  feature_info$t_statistic_direction <- sapply(1:num_features, function(i) ifelse(feature_info$t_statistic[i] > 0, "Up", "Down"))
  feature_info$pvalue <- sapply(1:num_features, function(i) t.test(cond_data[[2]][i, ], cond_data[[1]][i, ], var.equal = FALSE)$p.value)
  feature_info$qvalue <- p.adjust(feature_info$pvalue, "BH")
  feature_info$DEstatus <- sapply(1:num_features, function(i) ifelse(abs(feature_info$qvalue[i]) >= 0.05, FALSE, TRUE))

  return(feature_info)
}
