#'#'#'
#'
#'
#'#' Show function will display general information about the data present in the DNEAobject slots
#' @param object A DNEAobject
#' @export
#' @noRd
setMethod("show", "DNEAresults", function(object) {
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

expressionData.DNEAresults <- function(x, type = c("input", "normalized")){

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
setMethod("expressionData",signature(x = "DNEAresults"), expressionData.DNEAresults)

conditionLevels.DNEAresults <- function(x){
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
setMethod("conditionLevels", signature(x = "DNEAresults"), conditionLevels.DNEAresults)

conditions.DNEAresults <- function(x){

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
setMethod("conditions", signature(x = "DNEAresults"), conditions.DNEAresults)

#' @rdname conditions
#' @export
setMethod("conditions<-", signature(x = "DNEAresults"), function(x, value){

  x@metadata$samples$conditions <- x@metadata$samples[[value]]
  validObject(x)
  x
})

sampleNames.DNEAresults <- function(x){
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
setMethod("sampleNames", signature(x = "DNEAresults"), sampleNames.DNEAresults)

featureNames.DNEAresults <- function(x, original = TRUE){

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
setMethod("featureNames", signature(x = "DNEAresults"), featureNames.DNEAresults)

numFeatures.DNEAresults <- function(x){
  x@dataset_summary$num_features
}
#' numFeatures finds the total number of features in the dataset
#'
#' @param x A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return The number of features in the dataset
#'
#' @rdname numFeatures
#' @export
setMethod("numFeatures", signature(x = "DNEAresults"), numFeatures.DNEAresults)

numSamples.DNEAresults <- function(x){
  x@dataset_summary$num_samples
}
#' numSamples finds the total number of samples in the dataset
#'
#' @param x A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return The number of samples in the dataset
#' @rdname numSamples
#' @export
setMethod("numSamples", signature(x = "DNEAresults"), numSamples.DNEAresults)

optimizedLambda.DNEAresults <- function(x){
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
setMethod("optimizedLambda", signature(x = "DNEAresults"), optimizedLambda.DNEAresults)

#' @rdname optimizedLambda
#' @export
setReplaceMethod("optimizedLambda", signature(x = "DNEAresults"), function(x, value){

  x@hyperparameter$optimized_lambda <- value
  validObject(x)
  x
})

lambdas2Test.DNEAresults <- function(x){
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
setMethod("lambdas2Test", signature(x = "DNEAresults"), lambdas2Test.DNEAresults)

#' @rdname lambdas2Test
#' @export
setReplaceMethod("lambdas2Test", signature(x = "DNEAresults"), function(x, value){

  x@hyperparameter$tested_lambda_values <- value
  validObject(x)
  x
})

BICscores.DNEAresults <- function(x){
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
setMethod("BICscores", signature(x = "DNEAresults"), function(x) BICscores.DNEAresults)

#' @rdname BICscores
#' @keywords internal
setReplaceMethod("BICscores", signature(x = "DNEAresults"), function(x, value){

  x@hyperparameter$BIC_scores <- value
  validObject(x)
  x
})

selectionResults.DNEAresults <- function(x){
  x@stable_networks$selection_results
}
#' selectionResults returns an m x m matrix of the edge selection results from stabilitySelection()
#'
#' The function takes as input a pcorNetwork, DNEAobject, or DNEAobject_collapsed object and returns
#' an m x m matrix of selection results for every possible network edge calculated through
#' stabilitySelection()
#'
#' @param x A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return A pcorNetwork, DNEAobject, or DNEAobject_collapsed object after filling the
#' selection_results section of the stable_networks slot
#'
#' @rdname selectionResults
#' @keywords internal
setMethod("selectionResults", signature(x = "DNEAresults"), selectionResults.DNEAresults)

#'
#' @rdname selectionResults
#' @keywords internal
setReplaceMethod("selectionResults", signature(x = "DNEAresults"), function(x, value){

  x@stable_networks$selection_results <- value
  validObject(x)
  x
})

selectionProbabilities.DNEAresults <- function(x){
  x@stable_networks$selection_probabilities
}
#' selectionProbabilities returns an m x m matrix of the edge selection results from stabilitySelection()
#'
#' The function takes as input a pcorNetwork, DNEAobject, or DNEAobject_collapsed object and returns
#' an m x m matrix of selection Probabilities for every possible network edge calculated through
#' stabilitySelection()
#'
#' @param x A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return A pcorNetwork, DNEAobject, or DNEAobject_collapsed object after filling the
#' selection_probabilities section of the stable_networks slot
#'
#' @rdname selectionProbabilities
#' @keywords internal
setMethod("selectionProbabilities", signature(x = "DNEAresults"), selectionProbabilities.DNEAresults)

#'
#' @rdname selectionProbabilities
#' @keywords internal
setReplaceMethod("selectionProbabilities", signature(x = "DNEAresults"), function(x, value){

  x@stable_networks$selection_probabilities <- value
  validObject(x)
  x
})

edgeList.DNEAresults <- function(x){

  x@edge_list
}
#' edgeList returns the edge list created from DNEA
#'
#' The function takes as input a pcorNetwork, DNEAobject, or DNEAobject_collapsed object and returns
#'  the edge list created from the getNetworks() function
#'
#' @param x a pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return a data.frame corresponding to the edge list determined by DNEA
#'
#' @rdname edgeList
#' @export
setMethod("edgeList", signature(x = "DNEAresults"), edgeList.DNEAresults)

#' @rdname edgeList
#' @keywords internal
setReplaceMethod("edgeList", signature(x = "DNEAresults"), function(x, value){
  x@edge_list <- value
  validObject(x)
  x
})

nodeList.DNEAresults <- function(x){

  x@node_list
}
#' nodeList returns the node list created from DNEA
#'
#' The function takes as input a pcorNetwork, DNEAobject, or DNEAobject_collapsed object and returns
#'  the node list created from the getNetworks() function
#'
#' @param x a pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return a data.frame corresponding to the node list determined by DNEA
#'
#' @rdname nodeList
#' @export
setMethod("nodeList", signature(x = "DNEAresults"), nodeList.DNEAresults)

#' @rdname nodeList
#' @keywords internal
setReplaceMethod("nodeList", signature(x = "DNEAresults"), function(x, value){

  x@node_list <- value
  validObject(x)
  x
})
datasetSummary.DNEAresults <- function(x){

  message(paste0("Number of samples: ", numSamples(object)))
  message(paste0("Number of features: ", numFeatures(x)))
  message("Condition levels:")
  message(paste0("  1. ", conditionLevels(x)[1]))
  message(paste0("  2. ", conditionLevels(x)[2]))
  message("The diagnostic values are as follows:")
  print(as.matrix(x@dataset_summary$diagnostic_values))

  return(x@dataset_summary)
}
#' datasetSummary retrieves the dataset summary
#'
#' FOR INTERNAL USE ONLY - The function takes as input a pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' and returns the datasetSummary containing the number of features and samples, and the diagnostic values
#'
#' @param a pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return prints the values to console and returns a list of the aforementioned variables
#'
#' @rdname datasetSummary
#' @keywords internal
setMethod("datasetSummary", signature(x = "DNEAresults"), datasetSummary.DNEAresults)

#' @rdname datasetSummary
#' @keywords internal
setReplaceMethod("datasetSummary", signature(x = "DNEAresults"), function(x, value){

  x@dataset_summary <- value
  validObject(x)
  x
})

adjacencyMatrix.DNEAresults <- function(x, weighted = FALSE){

  if(weighted){

    x@adjacency_matrix$weighted_adjacency

  } else if(!weighted){

    x@adjacency_matrix$unweighted_adjacency
  }
}
#' adjacencyMatrix retrieves the weighted or unweighted adjacency matrix
#'
#' The function takes as input a pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' and returns the weighted or unweighted adjacency matrix determined via getNetworks().
#' **FOR INTERNAL USE ONLY - The matrix can also be inserted into the object via the setter functionality of the function.**
#'
#' @param a pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return The input object after filling the @@adjacency_matrix slot or a matrix corresponding to the adjacency matrix
#'
#' @rdname adjacencyMatrix
#' @keywords internal
setMethod("adjacencyMatrix", signature(x = "DNEAresults"), adjacencyMatrix.DNEAresults)

#' @rdname adjacencyMatrix
#' @keywords internal
setReplaceMethod("adjacencyMatrix", signature(x = "DNEAresults"), function(x, weighted = FALSE, value){

  if(weighted){

    x@adjacency_matrix$weighted_adjacency <- value

  } else if(!weighted){

    x@adjacency_matrix$unweighted_adjacency <- value
  }

  validObject(x)
  x
})

#' adjacencyGraph retrieves the adjacency graph for the case, control, or joint network
#'
#' The function takes as input a pcorNetwork, DNEAobject, DNEAobject_collapsed, or consensusClusteringResults object
#' and returns the adjacency graph made for the case, control, or joint network determined via runConsensusCluster().
#' **FOR INTERNAL USE ONLY - The graph can also be inserted into the object via the setter functionality of the function.**
#'
#' @param x A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @param graph the adjacency graph to retrieve. Values can be a provided condition name for either independent network,
#' or "joint_graph" for the total network graph
#' @return The input object after filling the @@adjacency_matrix slot or a graph corresponding to the specified
#' adjacency graph
#'
#' @rdname adjacencyGraph
#' @export
setMethod("adjacencyGraph", signature(x = "DNEAresults"), function(x, graph){

  x@consensus_clustering@adjacency_graphs$graph
})

#' @rdname adjacencyGraph
#' @keywords internal
setReplaceMethod("adjacencyGraph", signature(x = "DNEAresults"), function(x, graph, value){

  x@consensus_clustering@adjacency_graphs$graph <- value
  validObject(x)
  x
})

#' @rdname adjacencyGraph
#' @export
setMethod("adjacencyGraph", signature(x = "consensusClusteringResults"), function(x, graph){

  x@adjacency_graphs$graph
})

#' @rdname adjacencyGraph
#' @keywords internal
setReplaceMethod("adjacencyGraph", signature(x = "consensusClusteringResults"), function(x, graph, value){

  x@adjacency_graphs$graph <- value
  validObject(x)
  x
})

#' summary retrieves the summary results of consensus clustering
#'
#' The function takes as input a pcorNetwork, DNEAobject, DNEAobject_collapsed, or consensusClusteringResults object
#' and returns a summary of the results of consensus clustering determined via runConsensusCluster().
#'
#' @param x A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return A data.frame that corresponds to a summary of the results of consensus clustering
#'
#' @rdname summary
#' @export
setMethod("summary", signature(object = "consensusClusteringResults"), function(object){

  object@summary
})

#' @rdname summary
#' @export
setMethod("summary", signature(object = "DNEAresults"), function(object){

  object@consensus_clustering@summary
})

#' subnetworkMembership retrieves the subnetwork membership for each feature
#'
#' The function takes as input a pcorNetwork, DNEAobject, DNEAobject_collapsed, or consensusClusteringResults object
#' and returns the results of consensus clustering determined via runConsensusCluster().
#'
#' @param x A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return A data.frame that corresponds to the results of consensus clustering
#'
#' @rdname subnetworkMembership
#' @export
setMethod("subnetworkMembership", signature(x = "consensusClusteringResults"), function(x){

  x@subnetwork_membership
})

#' @rdname subnetworkMembership
#' @export
setMethod("subnetworkMembership", signature(x = "DNEAresults"), function(x){

  x@consensus_clustering@subnetwork_membership
})








