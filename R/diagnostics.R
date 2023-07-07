#' Calculate diagnostic criteria to determine stability of dataset
#'
#' This function takes as input a DNEAresults object and first splits the scaled data by condition. The minimum eigen
#' value and condition number are then calculated for the whole dataset as well as for each condition to determine
#' mathematic stability of the dataset and subsequent results from a GGM model. More information about interpretation can be
#' found in *Details*.
#'
#' @param mat A matrix of expression data
#' @param condition_values A vector of condition names present in the dataset
#' @param conditions A vector of condition labels corresponding to the samples in mat
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{createDNEAobject}}, \code{\link{reduceFeatures}}
#'
#' @details
#' Negative or zero eigenvalues in a dataset can represent instability in that portion of the matrix, respectively, thereby invalidating
#' parametric statistical methods and creating unreliable results. In this function, the minimum eigenvalue of the dataset
#' is calculated by first creating a pearson correlation matrix of the data. Negative eigenvalues may then occur for a
#' number of reasons, but one commonly seen cause is highly correlated features (in the positive and negative
#' direction). Regularization often takes care of this problem by arbitrarily selecting one of the variables in a highly
#' correlated group and removing the rest. Due to the *p >> n* problem common in -omics datasets, we have developed DNEA
#' to be very robust in removing redundant features via several regularization steps (**please see** \code{\link{BICtune}}
#' **and** \code{\link{stabilitySelection}}). As such, highly correlated features may already be handled downstream, however,
#' this may not be ideal for several reasons.
#'
#' In scenarios like this we recommend collapsing highly correlated features into a single group using
#' \code{\link{reduceFeatures}} - particularly if the dataset contains many highly-correlated features of a given class of molecules,
#' ie. many fatty acids or carnitines, respectively. Doing so gives the user more control over what variables are included
#' in the model. Without collapsing, the model regularization may result in one of the features within a class being included
#' and some or all of the remaining features being removed. By collapsing first, you retain the signal from all of the features
#' in the model and also have information pertaining to which features are highly correlated and as a result track each other.
#'
#' @returns A DNEA object containing an initialized node_list. Differential Expression is performed on
#'          the features if un-scaled data is provided. The min eigen value and condition number is also
#'          printed for the whole dataset as well as each condition.
#'
#' @importFrom stats t.test cor p.adjust
#' @keywords internal
#' @noRd
dataDiagnostics <- function(mat, condition_values, conditions) {

  ##set necessary parameters
  num_features <- ncol(mat)
  num_samples <- nrow(mat)

  ##split data by condition
  separated_conditions_scaled <- split_by_condition(dat = mat,
                                                    condition_levels = condition_values,
                                                    condition_by_sample = conditions)

  ##min eigen value and condition number for each of the datasets
  #control
  pearson <- cor(t(separated_conditions_scaled[[1]]))
  cond_number_1 <- kappa(pearson, exact=TRUE)
  min_eigen_1 <- min(unlist(eigen(pearson, symmetric=TRUE, only.values=TRUE)))

  #case
  pearson <- cor(t(separated_conditions_scaled[[2]]))
  cond_number_2 <- kappa(pearson, exact=TRUE)
  min_eigen_2 <- min(unlist(eigen(pearson, symmetric=TRUE, only.values=TRUE)))

  #whole dataset
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

  return(list(num_samples = num_samples, num_features = num_features,
              condition_levels = condition_values, diagnostic_values = diagnostic_values))
}




