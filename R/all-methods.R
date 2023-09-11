#'
#'
#'
#' Display general information about the data present in the DNEAresults object slots
#'
#' This function will display a summary of the information stored within a \code{DNEAresults} object
#'
#' @param object A \code{DNEAresults} object
#' @author Christopher Patsalis
#' @seealso \code{\link{createDNEAobject}}, \code{\link{reduceFeatures}}
#'
#' @returns A summary of the information stored in a \code{DNEAresults} object.
#' @examples
#' #import example data
#' data(dnw)
#'
#' dnw
#'
#' @keywords internal
#' @export
setMethod("show", "DNEAresults", function(object) {
  cat(is(object)[[1]], "\n",
      "  Project Name -  ", object@project_name, "\n",
      "  Samples -  ", paste0('There are ',length(sampleNames(object)), ' samples.'), "\n",
      "  Features -  ", paste0('There are ',length(featureNames(object)), ' Features.'), "\n",
      "  Conditions -  ", paste0("control: ", networkGroups(object)[[1]],
                                 " case: ", networkGroups(object)[[2]]), "\n",
      "  Optimized Lambda - ", paste0("lambda: ",optimizedLambda(object),
                                      " will be used in analysis."), "\n",
      "  Metabolic modules: ", paste0("there are ", length(unique(nodeList(object)[["membership"]])),
                                 " subnetworks."), "\n",
      "  netGSA results", paste0(sum(netGSAresults(object)$NetGSA_pFDR <= 0.05),
                                 " of the metabolic modules are differentially enriched across experimental condition"),
      sep = ""
  )
})

#' Display general information about the data present in the DNEAinputSummary slots
#'
#' This function will display the number of samples, number of features, and diagnostics values of the input dataset to a
#' \code{DNEAresults} object
#'
#' @param object A DNEAinputSummary object
#' @author Christopher Patsalis
#' @seealso \code{\link{createDNEAobject}}, \code{\link{reduceFeatures}}
#'
#' @returns A summary of the input data to \code{\link{createDNEAobject}}.
#' @keywords internal
#' @export
setMethod("show", "DNEAinputSummary", function(object) {
  cat(paste0(is(object)[[1]], "\n",
             "  Number of Samples  -  ", numSamples(object), "\n",
             "  Number of Features  -  ", numFeatures(object), "\n"),
      sep = "")
  print(as.matrix(diagnostics(object)))
})

#' Return the name of the current experiment
#'
#' This function returns the name of the DNEAdev experiment
#'
#' @param x A \code{DNEAresults} object
#' @returns The name of the DNEAdev experiment.
#' @examples
#' #import example data
#' data(dnw)
#'
#' projectName(dnw)
#' @author Christopher Patsalis
#' @seealso \code{\link{createDNEAobject}}
#'
#' @rdname projectName
#' @export
setMethod("projectName", signature(x= "DNEAresults"), function(x){

  x@project_name
})

expressionData.DNEAresults <- function(x, normalized = TRUE){

  ##pull data to return
  if(normalized){

    output <- x@assays$scaled_expression_data
  }else if(!normalized){

    output <- x@assays$expression_data
  }

  return(output)
}
#' Access expression data within a DNEAresults object,
#'
#' This function accesses the expression data stored in the @@assays slot of the \code{DNEAresults} object. The output is an
#' \emph{n x m} matrix with one row for each sample and one column for each feature in the data.
#'
#'
#' @param x A \code{DNEAresults} object
#' @param normalized A boolean indicating whether the normalized or original input data should be returned
#'
#' @author Christopher Patsalis
#' @seealso \code{\link{createDNEAobject}}, \code{\link{reduceFeatures}}

#' @returns The expression matrix specified by the user.
#' @examples
#' #import example data
#' data(dnw)
#'
#' expressionData(dnw, normalized = TRUE)
#'
#' @rdname expressionData
#' @export
setMethod("expressionData",signature(x = "DNEAresults"), expressionData.DNEAresults)

networkGroupIDs.DNEAresults <- function(x){
  x@metadata$network_group_IDs
}
#' Access and set the experimental group labels utilized by DNEAdev
#'
#' This function accesses the experimental group labels for each sample stored in the @@metadata slot of a \code{DNEAresults} object.
#'
#' @param x A \code{DNEAresults} object
#' @param value a character string name corresponding to a column name of the sample metadata data frame
#' @author Christopher Patsalis
#' @seealso \code{\link{includeMetadata}}
#' @returns A vector of the unique condition labels.
#' @examples
#' #import example data
#' data(dnw)
#'
#' networkGroupIDs(dnw)
#' @rdname networkGroupIDs
#' @export
setMethod("networkGroupIDs", signature(x = "DNEAresults"), networkGroupIDs.DNEAresults)

#' Retrieve the unique group values of the experimental condition
#'
#' This function takes in a \code{DNEAresults} object and returns the unique group values of the experimental condition in the dataset
#'
#' @param x A DNEAobject, or DNEAobject_collapsed object
#' @author Christopher Patsalis
#' @seealso \code{\link{networkGroupIDs}}, \code{\link{createDNEAobject}}
#' @returns A vector of the condition values.
#' @examples
#' #import example data
#' data(dnw)
#'
#' networkGroups(dnw)
#' @rdname networkGroups
#' @export
setMethod("networkGroups", signature(x = "DNEAresults"), function(x){

  x@metadata$network_groups
})

sampleNames.DNEAresults <- function(x){
  x@metadata[["samples"]]$samples
}

#' Retrieve the sample names from the metadata slot.
#'
#' This function accesses the sample names stored in the @@metadata slot of the \code{DNEAresults} object.
#'
#' @param x A \code{DNEAresults} object
#' @author Christopher Patsalis
#' @seealso \code{\link{createDNEAobject}}
#' @returns A character vector of sample names.
#' @examples
#' #import example data
#' data(dnw)
#'
#' sampleNames(dnw)
#' @rdname sampleNames
#' @export
setMethod("sampleNames", signature(x = "DNEAresults"), sampleNames.DNEAresults)

featureNames.DNEAresults <- function(x, original = FALSE){

  if(isTRUE(original)){
    x@metadata[["features"]]$feature_names
  }else if(isFALSE(original)){
    x@metadata[["features"]]$clean_feature_names
  }
}
#' Retrieve the feature names from the metadata slot.
#'
#' This function accesses the feature names stored in the @@metadata slot of the \code{DNEAresults} object.
#'
#' @param x A \code{DNEAresults} object
#' @param original "TRUE" returns the original feature names and "FALSE" returns the feature
#'  names that have been modified to avoid errors as a result of special characters.
#' @author Christopher Patsalis
#' @seealso \code{\link{createDNEAobject}}
#' @returns A character vector of feature names.
#' @examples
#' #import example data
#' data(dnw)
#'
#' featureNames(dnw, original = TRUE)
#' @rdname featureNames
#' @export
setMethod("featureNames", signature(x = "DNEAresults"), featureNames.DNEAresults)

numFeatures.DNEAresults <- function(x){
  x@dataset_summary@num_features
}
#' Retrieve the total number of features in the dataset
#'
#' This function prints to console the total number of features in the dataset
#'
#' @param x A \code{DNEAresults} object
#' @author Christopher Patsalis
#' @seealso \code{\link{createDNEAobject}}
#' @returns The number of features in the dataset.
#' @examples
#' #import example data
#' data(dnw)
#'
#' numFeatures(dnw)
#' @rdname numFeatures
#' @export
setMethod("numFeatures", signature(x = "DNEAresults"), numFeatures.DNEAresults)

#' @rdname numFeatures
#' @export
setMethod("numFeatures", signature(x = "DNEAinputSummary"), function(x){

  x@num_features
})
numSamples.DNEAresults <- function(x){
  x@dataset_summary@num_samples
}
#' Retrieves the total number of samples in the dataset
#'
#' This function prints to console the total number of samples in the dataset
#' @param x A \code{DNEAresults} object
#' @author Christopher Patsalis
#' @seealso \code{\link{createDNEAobject}}
#' @returns The number of samples in the dataset.
#' @examples
#' #import example data
#' data(dnw)
#'
#' numSamples(dnw)
#' @rdname numSamples
#' @export
setMethod("numSamples", signature(x = "DNEAresults"), numSamples.DNEAresults)

#' @rdname numSamples
#' @export
setMethod("numSamples", signature(x = "DNEAinputSummary"), function(x){

  x@num_samples
})
optimizedLambda.DNEAresults <- function(x){

  x@hyperparameter$optimized_lambda
}

#' Access the lambda value used in analysis
#'
#' The function takes as input a \code{DNEAresults} object and returns the hyperparameter (lambda) that is currently being
#' used for the analysis. The user may also provide a single-value numeric vector to change the lambda value for analysis
#'
#' @param x A \code{DNEAresults} object
#' @param value a single-value numeric vector corresponding to the lambda value to use in analysis
#' @author Christopher Patsalis
#' @seealso \code{\link{BICtune}}
#' @returns The optimized lambda hyperparameter.
#' @examples
#' #import example data
#' data(dnw)
#'
#' optimizedLambda(dnw)
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
#' Access the lambda values tested during hyperparameter optimization
#'
#' The function takes as input a \code{DNEAresults} object and returns the lambda values that were testing
#' during hyperparameter optimization performed via \code{\link{BICtune}}. \cr
#'
#'
#' @param x A pcorNetwork, DNEAobjct, or DNEAobject_collapsed object
#' @author Christopher Patsalis
#' @seealso \code{\link{BICtune}}
#' @returns The lambda values to evaluate in optimization.
#' @examples
#' #import example data
#' data(dnw)
#'
#' lambdas2Test(dnw)
#' @rdname lambdas2Test
#' @export
setMethod("lambdas2Test", signature(x = "DNEAresults"), lambdas2Test.DNEAresults)

#'
#' @rdname lambdas2Test
#' @keywords internal
setReplaceMethod("lambdas2Test", signature(x = "DNEAresults"), function(x, value){

  x@hyperparameter$tested_lambda_values <- value
  validObject(x)
  x
})

BICscores.DNEAresults <- function(x){

  x@hyperparameter$BIC_scores
}
#' Access the BIC scores for each lambda value evaluated
#'
#' The function takes as input a \code{DNEAresults} object and returns the BIC values for each lambda tested during
#' hyperparameter optimization performed via BICtune(). \cr
#'
#'
#' @param x A \code{DNEAresults} object
#' @author Christopher Patsalis
#' @seealso \code{\link{BICtune}}
#' @returns The optimized lambda hyperparameter.
#' @examples
#' #import example data
#' data(dnw)
#'
#' BICscores(dnw)
#' @rdname BICscores
#' @export
setMethod("BICscores", signature(x = "DNEAresults"), BICscores.DNEAresults)

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
#' Access and set the edge selection results from stabilitySelection()
#'
#' The function takes as input a \code{DNEAresults} object and returns an m x m matrix of selection results for every
#' possible network edge calculated via \code{\link{stabilitySelection}}. \cr
#'
#'
#' @param x A \code{DNEAresults} object
#' @author Christopher Patsalis
#' @seealso \code{\link{stabilitySelection}}, \code{\link{selectionProbabilities}}
#' @returns A \code{DNEAresults} object after filling the
#'         selection_results section of the stable_networks slot.
#' @examples
#' #import example data
#' data(dnw)
#'
#' selectionResults(dnw)
#' @rdname selectionResults
#' @export
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
#' Access and set the edge selection probabilities from stabilitySelection()
#'
#' The function takes as input a \code{DNEAresults} object and returns an m x m matrix of selection probabilities for every
#' possible network edge calculated via \code{\link{stabilitySelection}}. \cr
#'
#'
#' @param x A \code{DNEAresults} object
#' @author Christopher Patsalis
#' @seealso \code{\link{stabilitySelection}}, \code{\link{selectionResults}}
#' @returns A \code{DNEAresults} object after filling the
#'         selection_probabilities section of the stable_networks slot.
#' @examples
#' #import example data
#' data(dnw)
#'
#' selectionProbabilities(dnw)
#' @rdname selectionProbabilities
#' @export
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
#' Access the edge list
#'
#' The function takes as input a \code{DNEAresults} object and returns the edge list created from \code{\link{getNetworks}}. \cr
#'
#'
#' @param x a \code{DNEAresults} object
#' @author Christopher Patsalis
#' @seealso \code{\link{getNetworks}}, \code{\link{filterNetworks}}, \code{\link{getNetworkFiles}}
#' @returns A data frame corresponding to the edge list determined by DNEAdev.
#' @examples
#' #import example data
#' data(dnw)
#'
#' edgeList(dnw)
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

#' Access the node list
#'
#' The function takes as input a \code{DNEAresults} object and returns the node list created from \code{\link{createDNEAobject}}. \cr
#'
#'
#' @param x a \code{DNEAresults} object
#' @author Christopher Patsalis
#' @seealso \code{\link{createDNEAobject}}, \code{\link{clusterNet}}, \code{\link{getNetworkFiles}}
#' @returns A data frame corresponding to the node list determined by DNEAdev.
#' @examples
#' #import example data
#' data(dnw)
#'
#' nodeList(dnw)
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

diagnostics.DNEAresults <- function(x){

  x@dataset_summary@diagnostic_values
}
#' Retrieve the diagnostic values for the input expression data
#'
#' This function retrieves teh diagnostic values calculated for the input expression data to \code{\link{createDNEAobject}} \cr
#'
#'
#' @param a \code{DNEAresults} object or \code{DNEAinputSummary} object
#' @author Christopher Patsalis
#' @seealso \code{\link{createDNEAobject}}, \code{\link{reduceFeatures}}
#' @returns Returns the diagnostic values for the input expression data.
#' @examples
#' #import example data
#' data(dnw)
#'
#' diagnostics(dnw)
#' @rdname diagnostics
#' @export
setMethod("diagnostics", signature(x = "DNEAresults"), diagnostics.DNEAresults)

#' @rdname diagnostics
#' @keywords internal
setReplaceMethod("diagnostics", signature(x = "DNEAresults"), function(x, value){

  x@dataset_summary$diagnostic_values <- value
  validObject(x)
  x
})

#' @rdname diagnostics
setMethod("diagnostics", signature(x = "DNEAinputSummary"), function(x){

  x@diagnostic_values
})

#' @rdname diagnostics
#' @keywords internal
setReplaceMethod("diagnostics", signature(x = "DNEAinputSummary"), function(x, value){

  x@diagnostic_values <- value
  validObject(x)
  x
})

#' Access the dataset_summary slot of a DNEAresults object
#'
#' This function prints to console the number of samples, number of features, and diagnostic values of the input data to
#' \code{\link{createDNEAobject}} at initiation of the DNEAdev workflow. \cr
#'
#'
#' @param x A \code{DNEAresults} object
#' @author Christopher Patsalis
#' @seealso \code{\link{createDNEAobject}}, \code{\link{reduceFeatures}}
#' @returns The numbers of samples/features and diagnostic values of the input data calculated by
#'          \code{\link{createDNEAobject}}.
#' @examples
#' #import example data
#' data(dnw)
#'
#' datasetSummary(dnw)
#' @rdname datasetSummary
#' @export
setMethod("datasetSummary", signature(x = "DNEAresults"), function(x){

  x@dataset_summary
})

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
#' Retrieve the weighted or unweighted adjacency matrix
#'
#' The function takes as input a \code{DNEAresults} object and returns the weighted or unweighted adjacency matrix
#' determined via \code{\link{getNetworks}}. \cr
#'
#'
#' @param a \code{DNEAresults} object
#' @author Christopher Patsalis
#' @seealso \code{\link{getNetworks}}
#' @returns A matrix corresponding to the adjacency matrix specified.
#' @examples
#' #import example data
#' data(dnw)
#'
#' adjacencyMatrix(dnw, weighted = TRUE)
#' @rdname adjacencyMatrix
#' @export
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

#' Retrieve the adjacency graph for the case, control, or joint network
#'
#' The function  returns the adjacency graph made for the case, control, or joint network determined via
#' \code{\link{clusterNet}}. \cr
#'
#'
#' @param x A \code{DNEAresults} or \code{consensusClusteringResults} object
#' @param graph A character string indicating which of the adjacency graphs to return. Values can be "joint_graph" for the
#' whole graph object, or one of the group values returned by \code{\link{networkGroups}}
#' @author Christopher Patsalis
#' @seealso \code{\link{clusterNet}}
#' @returns An \code{\link{igraph}} graph object corresponding to the specified adjacency graph.
#' @examples
#' #import example data
#' data(dnw)
#'
#' adjacencyGraph(dnw, graph = "DM:case")
#' @rdname adjacencyGraph
#' @export
setMethod("adjacencyGraph", signature(x = "DNEAresults"), function(x, graph){

  x@consensus_clustering@adjacency_graphs[[graph]]
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

#' Retrieve a summary of a consensusClusteringResults object
#'
#' The function takes as input a consensusClusteringResults object and returns a summary of the results of
#' consensus clustering determined via \code{\link{clusterNet}}.
#'
#' @param object A consensusClusteringResults
#' @author Christopher Patsalis
#' @seealso \code{\link{clusterNet}}
#' @returns A data frame that corresponds to a summary of the results of consensus clustering.
#' @rdname summary
#' @export
setMethod("summary", signature(object = "consensusClusteringResults"), function(object){

  object@summary
})

#' Retrieves the summary results of consensus clustering
#'
#' The function takes as input a \code{DNEAresults} object and returns a summary
#' of the results of consensus clustering determined via \code{\link{clusterNet}}.
#'
#' @param x A \code{DNEAresults} object
#' @author Christopher Patsalis
#' @seealso \code{\link{clusterNet}}
#' @returns A data frame summary of the consensus clustering results from DNEAdev.
#' @examples
#' #import example data
#' data(dnw)
#'
#' CCsummary(dnw)
#' @rdname CCsummary
#' @export
setMethod("CCsummary", signature(x = "DNEAresults"), function(x){

  summary(x@consensus_clustering)
})

#' @rdname CCsummary
#' @keywords internal
setReplaceMethod("CCsummary", signature(x = "DNEAresults"), function(x, value){

  x@consensus_clustering@summary <- value
  validObject(x)
  x
})
#' Retrieve the subnetwork membership for each feature
#'
#' The function takes as input a \code{DNEAresults} or \code{consensusClusteringResults} object
#' and returns the results of consensus clustering determined via \code{\link{clusterNet}}.
#'
#' @param x A \code{DNEAresults} object
#' @author Christopher Patsalis
#' @seealso \code{\link{clusterNet}}
#' @returns A data frame that corresponds to the results of consensus clustering.
#' @examples
#' #import example data
#' data(dnw)
#'
#' subnetworkMembership(dnw)
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

#' @rdname subnetworkMembership
#' @keywords internal
setReplaceMethod("subnetworkMembership", signature(x = "DNEAresults"), function(x, value){

  x@consensus_clustering@subnetwork_membership <- value
  validObject(x)
  x
})
#' Access the netGSA slot of a DNEAresults object
#'
#' The function takes as input a \code{DNEAresults} object and returns the netGSA results in the netGSA slot. \cr
#'
#' @param x A \code{DNEAresults} object
#' @author Christopher Patsalis
#' @seealso \code{\link{runNetGSA}}, \code{\link[netgsa:NetGSA]{netgsa::NetGSA()}}
#' @returns A data frame of the results from netGSA.
#' @examples
#' #import example data
#' data(dnw)
#'
#' netGSAresults(dnw)
#' @rdname netGSAresults
#' @export
setMethod("netGSAresults", signature(x = "DNEAresults"), function(x){

  x@netGSA
})

#' @rdname netGSAresults
#' @keywords internal
setReplaceMethod("netGSAresults", signature(x = "DNEAresults"), function(x, value){

  x@netGSA <- value
  validObject(x)
  x
})
