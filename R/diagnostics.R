#' Datadiagnostics will run Differential Expression on Nodes and check data stability prior to analysis
#'
#' This function first takes in a DNEAobject and runs diagnostics on the scaled expression data in the
#' scaled_expression_data value of the assays tab. It will then calculate the minimum eigen value and
#' condition number for the whole dataset, as well as for each condition and inform the user if they can
#' proceed with the analysis or collapse the features using the reduceFeatures() function prior.
#'
#' @param object A DNEAobject
#'
#' @returns A DNEA object containing an initialized node_list. Differential Expression is performed on
#'          the features if un-scaled data is provided. The min eigen value and condition number is also
#'          printed for the whole dataset as well as each condition.
#'
#' @importFrom stats t.test cor p.adjust
#'@noRd
dataDiagnostics <- function(object) {

  ############################################
  #**Set-up output and initialize variables**#
  ############################################

  num_features <- length(featureNames(object))
  num_samples <- length(sampleNames(object))

  #Set up output structure
  feature_info <- data.frame('Features' = featureNames(object),
                             'clean_Feature_Names' = cleanFeatureNames(object))

  #create list of conditions in data
  condition_levels <- unlist(levels(condition(object)))
  names(condition_levels) <- condition_levels



  ######################################
  #**Differential Expression Analysis**#
  ######################################

  #DE can only be done on un-scaled data because scaling each condition separately results in a mean
  #of zero for every feature and no difference in calculations.
  if(is.null(expressionData(object)) == FALSE){

    #split data by condition to perform DE
    cond_data <- split_by_condition(dat = expressionData(object),
                                    condition_levels = condition_levels,
                                    condition_by_sample = condition(object))

    ##fill in DE in node_lsit
    assign(paste0('feature_info$foldchange-',levels(condition(object))[2], '/', levels(condition(object))[1]), rowMeans(cond_data[[2]]) - rowMeans(cond_data[[1]]))
    feature_info$fcdirection <- sapply(1:num_features, function(i) ifelse(get(paste0('feature_info$foldchange-',levels(condition(object))[2], '/', levels(condition(object))[1]))[i] > 0, "Up", "Down"))
    feature_info$t_statistic <- sapply(1:num_features, function(i) t.test(cond_data[[2]][i, ], cond_data[[1]][i, ], var.equal = FALSE)$statistic)
    feature_info$t_statistic_direction <- sapply(1:num_features, function(i) ifelse(feature_info$t_statistic[i] > 0, "Up", "Down"))
    feature_info$pvalue <- sapply(1:num_features, function(i) t.test(cond_data[[2]][i, ], cond_data[[1]][i, ], var.equal = FALSE)$p.value)
    feature_info$qvalue <- p.adjust(feature_info$pvalue, "BH")
    feature_info$DEstatus <- sapply(1:num_features, function(i) ifelse(abs(feature_info$qvalue[i]) >= 0.05, FALSE, TRUE))

  } else{

    cat("No un-scaled data - Differential Expression on features (ie. nodes) could not be performed. \n")
    cat("Please see the vignette for more information on why this has ocurred. \n\n")
  }

  ######################################
  #**Feature Diagnostics**#
  ######################################

  separated_conditions_scaled <- split_by_condition(dat = scaledExpressionData(object),
                                                    condition_levels = condition_levels,
                                                    condition_by_sample = condition(object))
  #Calculate important values
  pearson <- cor(t(separated_conditions_scaled[[1]]))
  cond_number_1 <- kappa(pearson, exact=TRUE)
  min_eigen_1 <- min(unlist(eigen(pearson, symmetric=TRUE, only.values=TRUE)))
  pearson <- cor(t(separated_conditions_scaled[[2]]))
  cond_number_2 <- kappa(pearson, exact=TRUE)
  min_eigen_2 <- min(unlist(eigen(pearson, symmetric=TRUE, only.values=TRUE)))
  pearson <- cor(t(cbind(separated_conditions_scaled[[1]], separated_conditions_scaled[[2]])))
  cond_number_3 <- kappa(pearson, exact=TRUE)
  min_eigen_3 <- min(unlist(eigen(pearson, symmetric=TRUE, only.values=TRUE)))

  #create diagnostics table
  diagnostic_values <- data.frame(row.names = c("all_data",
                                                levels(condition(object))[[1]],
                                                levels(condition(object))[[2]]))
  diagnostic_values$min_eigen <- list(min_eigen_3,min_eigen_1, min_eigen_2)
  diagnostic_values$condition_num <- list(cond_number_3, cond_number_1, cond_number_2)

  #print values and suggestions
  cat('Diagnostic values are as follows: \n')
  print(as.matrix(diagnostic_values))
  cat('\n')
  if(any(diagnostic_values$min_eigen <=1e-5)){
    warning('One or more conditions look unstable. You should collapse features before continuing the analysis!\n')
  } else{
    cat('Diagnostic tests complete. You can proceed with the analysis!\n')
  }

  return(list(list(num_samples = num_samples, num_features = num_features, diagnostic_values = diagnostic_values, condition_levels = condition_levels),feature_info))
}




