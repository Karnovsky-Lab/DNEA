#'
#'
#'
#'#' Show function will display general information about the data present in the DNEAobject slots
#' @param object A DNEAobject
#' @export
#' @noRd
setMethod("show", "pcorNetwork", function(object) {
  cat(is(object)[[1]], "\n",
      "  Project Name -  ", object@project_name, "\n",
      "  Un-scaled data - ", class(object@assays[["expression_data"]])[[1]], "\n",
      "  Scaled data -  ", class(object@assays[["scaled_expression_data"]])[[1]], "\n",
      "  Samples -  ", paste0('There are ',length(object@metadata[["samples"]]), ' samples.'), "\n",
      "  Features -  ", paste0('There are ',length(object@metadata[["features"]]), ' Features.'), "\n",
      "  Conditions -  ", paste0("control: ", object@dataset_summary[["condition_levels"]][[1]],
                                 " case: ", object@dataset_summary[["condition_levels"]][[1]]), "\n",
      "  Optimized Lambda - ", paste0("lambda: ",object@hyperparameter[["optimizedLambda"]],
                                      " will be used in analysis."), "\n",
      "  Sub-clusters: ", paste0("there are ", length(unique(object@node_list[["membership"]])),
                                 " sub-clusters."),
      sep = ""
  )
})

expressionData.pcorNetwork <- function(x, type = c("input", "normalized")){

  ##set data to be pulled
  type = match.arg(type)

  ##pull data to return
  if(type == "input"){
    output <- x@assays$expression_data
  }else if(type == "normalized"){
    output <- x@assays$scaled_expression_data
  }

  return(output)
}
#' expressionData is an accessor function for the data in the "assays" slot of a pcorNetwork,
#'  DNEAobject, or DNEAobject_collapsed
#'
#' The assays slot hold the expression data in a n x m matrix, with one row for each sample, and
#'  one column for each feature.
#'
#'
#' @param x A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @param type "input" will return the original data, and "normalized" will return the scaled data
#' @return The indicated expression expression matrix
#'
#' @rdname expressionData
#' @export
setMethod("expressionData",signature(x = "pcorNetwork"), expressionData.pcorNetwork)

conditionLevels.DNEAobject <- function(x){
  x@dataset_summary$condition_levels
}
#' conditionLevels returns the group information for samples in a DNEAobject,
#'  or DNEAobject_collapsed
#'
#' This function takes in a DNEAobject or DNEAobject_collapsed and return the unique condition
#'  labels located in the dataset_summary slot of the object
#'
#'
#' @param x A DNEAobject or DNEAobject_collapsed
#' @return A vector of the unique condition labels
#'
#' @rdname conditionLevels
#' @export
setMethod("conditionLevels", signature(x = "DNEAobject"), conditionLevels.DNEAobject)

conditions.DNEAobject <- function(x){

  x@metadata$samples$conditions
}
#' condition retrieves the condition values for each sample from the metadata slot.
#'
#' This function takes in a DNEAobject, or DNEAobject_collapsed object and
#'  return the condition values located in the metatdata slot of the object
#'
#' @param x A DNEAobject, or DNEAobject_collapsed object
#' @return A vector of the condition values
#'
#' @rdname conditions
#' @export
setMethod("conditions", signature(x = "DNEAobject"), function(x) conditions.DNEAobject)

#' @rdname conditions
#' @export
setMethod("conditions<-", signature(x = "DNEAobject"), function(x, value){

  x@metadata$samples$conditions <- x@metadata$samples[[value]]
  validObject(x)
  x
})

sampleNames.pcorNetwork <- function(x){
  x@metadata[["samples"]]$samples
}
#' sampleNames retrieves the sample names from the metadata slot.
#'
#' This function takes in a pcorNetwork, DNEAobject, or DNEAobject_collapsed object and
#'  return the sample names located in the metadata slot of the object.
#'
#' @param x A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return A vector of sample names
#'
#' @rdname sampleNames
#' @export
setMethod("sampleNames", signature(x = "pcorNetwork"), sampleNames.pcorNetwork)

featureNames.pcorNetwork <- function(x, original = TRUE){

  if(isTRUE(original)){
    x@metadata[["features"]]$feature_names
  }else if(isFALSE(original)){
    x@metadata[["features"]]$clean_feature_names
  }
}
#' featureNames retrieves the feature names from the metadata slot.
#'
#' This function takes in a pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#'  and return the specified feature names located in the metadata slot of the object.
#'
#' @param x A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @param original "TRUE" returns the original feature names and "FALSE" returns the feature
#'  names that have been modified to avoid errors as a result of special characters.
#'
#' @return A vector of feature names
#' @rdname featureNames
#' @export
setMethod("featureNames", signature(x = "pcorNetwork"), featureNames.pcorNetwork)

numFeatures.pcorNetwork <- function(x){
  x@dataset_summary$num_features
}
#' numFeatures finds the total number of features in the dataset
#'
#' @param x A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return The number of features in the dataset
#'
#' @rdname numFeatures
#' @export
setMethod("numFeatures", signature(x = "pcorNetwork"), numFeatures.pcorNetwork)

numSamples.pcorNetwork <- function(x){
  x@dataset_summary$num_samples
}
#' numSamples finds the total number of samples in the dataset
#'
#' @param x A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return The number of samples in the dataset
#' @rdname numSamples
#' @export
setMethod("numSamples", signature(x = "pcorNetwork"), function(x) numSamples.pcorNetwork)

optimizedLambda.pcorNetwork <- function(x){
  x@hyperparameter$optimized_lambda
}
#' optimizedLambda returns the lambda value used in analysis.
#'
#' The function takes as input a pcorNetwork, DNEAobject, or DNEAobject_collapsed object and returns
#'  the hyperparameter (lambda) that is currently being used for the analysis
#'
#' @param x A pcorNetwork, DNEAobjct, or DNEAobject_collapsed object
#' @return The optimized lambda hyperparameter
#'
#' @rdname optimizedLambda
#' @export
setMethod("optimizedLambda", signature(x = "pcorNetwork"), optimizedLambda.pcorNetwork)

#' @rdname optimizedLambda
#' @export
setReplaceMethod("optimizedLambda", signature(x = "pcorNetwork"), function(x, value){

  x@hyperparameter$optimized_lambda <- value
  validObject(x)
  x
})

lambdas2Test.pcorNetwork <- function(x){
  x@hyperparameter$tested_lambda_values
}
#' lambdas2Test returns the lambda values tested during hyperparameter optimization
#'
#' The function takes as input a pcorNetwork, DNEAobject, or DNEAobject_collapsed object and returns
#'  the lambda values that were testing during hyperparameter optimization performed via BICtune()
#'
#' @param x A pcorNetwork, DNEAobjct, or DNEAobject_collapsed object
#' @return The lambda values to evaluate in optimization
#'
#' @rdname lambdas2Test
#' @export
setMethod("lambdas2Test", signature(x = "pcorNetwork"), lambdas2Test.pcorNetwork)

#' @rdname lambdas2Test
#' @export
setReplaceMethod("lambdas2Test", signature(x = "pcorNetwork"), function(x, value){

  x@hyperparameter$tested_lambda_values <- value
  validObject(x)
  x
})

BICscores.pcorNetwork <- function(x){
  x@hyperparameter$BIC_scores
}
#' BICscores returns the BIC scores for each lambda value evaluated
#'
#' The function takes as input a pcorNetwork, DNEAobject, or DNEAobject_collapsed object and returns
#'  the BIC values for each lambda tested during hyperparameter optimization performed via BICtune().
#'
#' @param x A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return The optimized lambda hyperparameter
#'
#' @rdname BICscores
#' @export
setMethod("BICscores", signature(x = "pcorNetwork"), function(x) BICscores.pcorNetwork)






