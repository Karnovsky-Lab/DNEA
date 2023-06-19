#' Datadiagnostics will run Differential Expression on Nodes and check data stability prior to analysis
#'
#' This function first takes in a DNEAobject and runs diagnostics on the scaled expression data in the
#' scaled_expression_data value of the assays tab. It will then calculate the minimum eigen value and
#' condition number for the whole dataset, as well as for each condition and inform the user if they can
#' proceed with the analysis or collapse the features using the reduceFeatures() function prior.
#'
#' @param mat A matrix of expression data
#' @param condition_values A vector of condition names present in the dataset
#' @param conditions A vector of condition labels corresponding to the samples in mat
#' @returns A DNEA object containing an initialized node_list. Differential Expression is performed on
#'          the features if un-scaled data is provided. The min eigen value and condition number is also
#'          printed for the whole dataset as well as each condition.
#'
#' @importFrom stats t.test cor p.adjust
#' @keywords internal
dataDiagnostics <- function(mat, condition_values, conditions) {

  ##set necessary parameters
  num_features <- ncol(mat)
  num_samples <- nrow(mat)

  ##split data by condition
  separated_conditions_scaled <- split_by_condition(dat = mat,
                                                    condition_levels = condition_values,
                                                    condition_by_sample = conditions)

  ##min eigen value and condition number for each of the datasets
  pearson <- cor(t(separated_conditions_scaled[[1]]))
  cond_number_1 <- kappa(pearson, exact=TRUE)
  min_eigen_1 <- min(unlist(eigen(pearson, symmetric=TRUE, only.values=TRUE)))
  pearson <- cor(t(separated_conditions_scaled[[2]]))
  cond_number_2 <- kappa(pearson, exact=TRUE)
  min_eigen_2 <- min(unlist(eigen(pearson, symmetric=TRUE, only.values=TRUE)))
  pearson <- cor(t(cbind(separated_conditions_scaled[[1]], separated_conditions_scaled[[2]])))
  cond_number_3 <- kappa(pearson, exact=TRUE)
  min_eigen_3 <- min(unlist(eigen(pearson, symmetric=TRUE, only.values=TRUE)))

  ##create diagnostics table
  diagnostic_values <- data.frame(row.names = c("all_data",
                                                condition_values[[1]],
                                                condition_values[[2]]))
  diagnostic_values$min_eigen <- list(min_eigen_3,min_eigen_1, min_eigen_2)
  diagnostic_values$condition_num <- list(cond_number_3, cond_number_1, cond_number_2)

  #print values and suggestions
  cat('Diagnostic values are as follows: \n')
  print(as.matrix(diagnostic_values))
  cat('\n')

  if(any(diagnostic_values$min_eigen <= 1e-5)){
    warning('One or more conditions look unstable. You should collapse features before continuing the analysis!\n')
  } else{
    cat('Diagnostic tests complete. You can proceed with the analysis!\n')
  }

  return(list(num_samples = num_samples, num_features = num_features, diagnostic_values = diagnostic_values, condition_levels = condition_values))
}




