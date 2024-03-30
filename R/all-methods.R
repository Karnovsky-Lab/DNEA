#' Retrieve information about a DNEAobj object
#'
#' @describeIn DNEAobj-class
#' This function will display a summary of the information stored within
#' a \code{\link{DNEAobj}} object
#'
#'
#' @param object A \code{\link{DNEAobj}} object
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}},
#' \code{\link{aggregateFeatures}}
#'
#' @returns A summary of the information stored in a
#' \code{\link{DNEAobj}} object.
#' @examples
#' #import example data
#' data(dnw)
#'
#' dnw
#' @export
setMethod("show", "DNEAobj", function(object) {
  cat(is(object)[[1]],
      "\n  Project Name -  ", object@project_name,
      "\n  Samples -  ", 'There are ',length(sampleNames(object)),
      ' samples.',
      "\n  Features -  ", 'There are ',length(featureNames(object)),
      ' Features.',
      "\n  Conditions -  ", "control: ", networkGroups(object)[[1]],
                            " case: ", networkGroups(object)[[2]],
      "\n  Optimized Lambda - ", ifelse(is.null(optimizedLambda(object)),
                                        "Parameter tuning not performed",
                                        paste0("lambda: ",
                                               optimizedLambda(object),
                                               " will be used in analysis.")),
      "\n  Metabolic modules: ", ifelse(is.null(nodeList(object)[["membership"]]),
                                        "Consensus Clustering not performed",
                                        paste0("there are ",
                                               length(unique(nodeList(object)[["membership"]])),
                                               " subnetworks.")),
      "\n  netGSA results: ", ifelse(nrow(netGSAresults(object)) == 0,
                                     "Enrichment analysis not performed",
                                     paste0(sum(netGSAresults(object)$NetGSA_pFDR <= 0.05),
                                            " enriched modules\n")),
      sep=""
  )
})

#' Retrieve information about a DNEAinputSummary object
#'
#' @describeIn DNEAinputSummary-class
#' This function will display the number of samples, number of features,
#' and diagnostics values of the input dataset to a
#' \code{\link{DNEAobj}} object.
#'
#' @param object A DNEAinputSummary object
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}},
#' \code{\link{aggregateFeatures}}
#'
#' @returns A summary of the input data to \code{\link{createDNEAobject}}.
#' @export
setMethod("show", "DNEAinputSummary", function(object){
  cat(is(object)[[1]],
             "\n  Number of Samples  -  ", numSamples(object),
             "\n  Number of Features  -  ", numFeatures(object),"\n",
      sep="")
  print(as.matrix(diagnostics(object)))
})

projectName.DNEAobj <- function(x){

  x@project_name
}
#' Return the name of the current experiment
#'
#' This function returns the name of the DNEA experiment
#'
#' @param x A \code{\link{DNEAobj}} object
#' @returns The name of the DNEA experiment.
#' @examples
#' #import example data
#' data(dnw)
#'
#' projectName(dnw)
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}}
#' @rdname projectName-methods
#' @aliases projectName
#' @export
setMethod("projectName", signature(x= "DNEAobj"),
          projectName.DNEAobj)

expressionData.DNEAobj <- function(x, assay=c("input_data",
                                              "log_input_data",
                                              "scaled_expression_data")){

  assay <- match.arg(assay)
  output <- x@assays[[assay]]
  return(output)
}
#' Access expression data within a DNEAobj object,
#'
#' This function accesses the expression data stored in the @@assays
#' slot of the \code{\link{DNEAobj}} object. The output is an \emph{n x m}
#' matrix with one row for each sample and one column for each feature
#' in the data.
#'
#'
#' @param x A \code{\link{DNEAobj}} object
#' @param assay A character string corresponding to the data to retrieve:
#' "input_data" retrieves the data as it was input, "log_input_data"
#' retrieves the input_data after log transforming, and
#' "scaled_expression_data" retrieves the log-scaled data matrices.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}},\code{\link{aggregateFeatures}}

#' @returns The expression matrix specified by the user.
#' @examples
#' #import example data
#' data(dnw)
#'
#' expressionData(x=dnw, assay="input_data")
#' @rdname expressionData-methods
#' @aliases expressionData
#' @export
setMethod("expressionData",signature(x="DNEAobj"),
          expressionData.DNEAobj)

assays.DNEAobj <- function(x){

  x@assays
}
#' @rdname assays-methods
#' @aliases assays
#' @keywords internal
#' @noRd
setMethod("assays", signature(x="DNEAobj"),
          assays.DNEAobj)

networkGroupIDs.DNEAobj <- function(x){

  x@metadata$network_group_IDs
}
#' Access and set the experimental group labels
#'
#' This function accesses the experimental group labels for each sample
#' stored in the @@metadata slot of a \code{\link{DNEAobj}} object.
#'
#' @param x A \code{\link{DNEAobj}} object
#' @param value a character string name corresponding to a column
#' name of the sample metadata data frame
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{includeMetadata}}
#' @returns A vector of the unique condition labels.
#' @examples
#' #import example data
#' data(dnw)
#'
#' networkGroupIDs(dnw)
#' @rdname networkGroupIDs-methods
#' @aliases networkGroupIDs
#' @export
setMethod("networkGroupIDs", signature(x="DNEAobj"),
          networkGroupIDs.DNEAobj)

networkGroups.DNEAobj <- function(x){

  x@metadata$network_groups
}
#' Retrieve the unique group values of the experimental condition
#'
#' This function takes in a \code{\link{DNEAobj}} object and returns the
#' unique group values of the experimental condition in the dataset
#'
#' @param x A DNEAobject, or DNEAobject_collapsed object
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{networkGroupIDs}},
#' \code{\link{createDNEAobject}}
#' @returns A vector of the condition values.
#' @examples
#' #import example data
#' data(dnw)
#'
#' networkGroups(dnw)
#' @rdname networkGroups-methods
#' @aliases networkGroups
#' @export
setMethod("networkGroups", signature(x="DNEAobj"),
          networkGroups.DNEAobj)

metaData.DNEAobj <- function(x, type=c("samples", "features")){

  type <- match.arg(type)
  if(type == "samples"){
    x@metadata[["samples"]]
  }else if(type == "features"){
    x@metadata[["features"]]
  }
}
#' Retrieve metadata stored in DNEAobj
#'
#' This function retrieves the specified metadata data frame stored
#' in the @@metadata slot of the \code{\link{DNEAobj}} object.
#'
#' @param x A \code{\link{DNEAobj}} object.
#' @param type A character string indicating the type of metadata
#' to access. Can be "sample" or "feature".
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}}, \code{\link{includeMetadata}}
#' @returns A data frame of the indicated metadata
#' @examples
#' #import example data
#' data(dnw)
#'
#' metaData(dnw, type = "sample")
#' @seealso \code{\link{includeMetadata}}
#' @rdname metaData-methods
#' @aliases metaData
#' @export
setMethod("metaData", signature(x="DNEAobj"),
          metaData.DNEAobj)

metaDataReplace.DNEAobj <- function(x, type, value){

  type <- match.arg(type)
  if(type == "samples"){
    x@metadata[["samples"]] <- value
  }else if(type == "features"){
    x@metadata[["features"]] <- value
  }
}
#' @keywords internal
#' @export
setReplaceMethod("metaData", signature(x="DNEAobj"),
                 metaDataReplace.DNEAobj)

sampleNames.DNEAobj <- function(x){

  x@metadata[["samples"]]$samples
}
#' Retrieve the sample names from the metadata slot.
#'
#' This function accesses the sample names stored in the @@metadata
#' slot of the \code{\link{DNEAobj}} object.
#'
#' @param x A \code{\link{DNEAobj}} object
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}}
#' @returns A character vector of sample names.
#' @examples
#' #import example data
#' data(dnw)
#'
#' sampleNames(dnw)
#' @rdname sampleNames-methods
#' @aliases sampleNames
#' @export
setMethod("sampleNames", signature(x="DNEAobj"),
          sampleNames.DNEAobj)

featureNames.DNEAobj <- function(x, original=FALSE){

  if(isTRUE(original)){
    x@metadata[["features"]]$feature_names
  }else if(isFALSE(original)){
    x@metadata[["features"]]$clean_feature_names
  }
}
#' Retrieve the feature names from the metadata slot.
#'
#' This function accesses the feature names stored in the @@metadata
#' slot of the \code{\link{DNEAobj}} object.
#'
#' @param x A \code{\link{DNEAobj}} object
#' @param original "TRUE" returns the original feature names and "FALSE"
#' returns the feature names that have been modified to avoid errors as
#' a result of special characters.
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}}
#' @returns A character vector of feature names.
#' @examples
#' #import example data
#' data(dnw)
#'
#' featureNames(dnw, original=TRUE)
#' @rdname featureNames-methods
#' @aliases featureNames
#' @export
setMethod("featureNames", signature(x="DNEAobj"),
          featureNames.DNEAobj)

numFeatures.DNEAobj <- function(x){
  x@dataset_summary@num_features
}
#' Retrieve the total number of features in the dataset
#'
#' This function prints to console the total number of features
#' in the dataset
#'
#' @param x A \code{\link{DNEAobj}} object
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}}
#' @returns The number of features in the dataset.
#' @examples
#' #import example data
#' data(dnw)
#'
#' numFeatures(dnw)
#' @rdname numFeatures-methods
#' @aliases numFeatures
#' @export
setMethod("numFeatures", signature(x="DNEAobj"),
          numFeatures.DNEAobj)

numFeatures.DNEAinputSummary <- function(x){

  x@num_features
}
#' @rdname numFeatures-methods
#' @export
setMethod("numFeatures", signature(x="DNEAinputSummary"),
          numFeatures.DNEAinputSummary)

numSamples.DNEAobj <- function(x){

  x@dataset_summary@num_samples
}
#' Retrieves the total number of samples in the dataset
#'
#' This function prints to console the total number of
#' samples in the dataset
#' @param x A \code{\link{DNEAobj}} object
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}}
#' @returns The number of samples in the dataset.
#' @examples
#' #import example data
#' data(dnw)
#'
#' numSamples(dnw)
#' @rdname numSamples-methods
#' @aliases numSamples
#' @export
setMethod("numSamples", signature(x="DNEAobj"),
          numSamples.DNEAobj)

numSamples.DNEAinputSummary <- function(x){

  x@num_samples
}
#' @rdname numSamples-methods
#' @aliases numSamples
#' @export
setMethod("numSamples", signature(x="DNEAinputSummary"),
          numSamples.DNEAinputSummary)

optimizedLambda.DNEAobj <- function(x){

  x@hyperparameter$optimized_lambda
}
#' Access the lambda value used in analysis
#'
#' The function takes as input a \code{\link{DNEAobj}} object and returns
#' the hyperparameter (lambda) that is currently being used for the analysis.
#' The user may also provide a single-value numeric vector to change the
#' lambda value for analysis
#'
#' @param x A \code{\link{DNEAobj}} object
#' @param value a single-value numeric vector corresponding to the lambda
#' value to use in analysis
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{BICtune}}
#' @returns The optimized lambda hyperparameter.
#' @examples
#' #import example data
#' data(dnw)
#'
#' optimizedLambda(dnw)
#' @rdname optimizedLambda-methods
#' @aliases optimizedLambda
#' @export
setMethod("optimizedLambda", signature(x="DNEAobj"),
          optimizedLambda.DNEAobj)

optimizedLambdaReplace.DNEAobj <- function(x, value){

  x@hyperparameter$optimized_lambda <- value
  validObject(x)
  x
}
#' @rdname optimizedLambda-methods
#' @aliases optimizedLambda
#' @export
setReplaceMethod("optimizedLambda", signature(x="DNEAobj"),
                 optimizedLambdaReplace.DNEAobj)

lambdas2Test.DNEAobj <- function(x){

  x@hyperparameter$tested_lambda_values
}
#' Access the lambda values tested during hyperparameter optimization
#'
#' The function takes as input a \code{\link{DNEAobj}} object and returns
#' the lambda values that were testing during hyperparameter optimization
#' performed via \code{\link{BICtune}}.
#'
#'
#' @param x A \code{\link{DNEAobj}} or \code{\link{collapsed_DNEAobj}} object
#' @param value A list or numeric vector
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{BICtune}}
#' @returns The lambda values to evaluate in optimization.
#' @examples
#' #import example data
#' data(dnw)
#'
#' lambdas2Test(dnw)
#' @rdname lambdas2Test-methods
#' @aliases lambdas2Test
#' @export
setMethod("lambdas2Test", signature(x="DNEAobj"),
          lambdas2Test.DNEAobj)

lambdas2TestReplace.DNEAobj <- function(x, value){

  x@hyperparameter$tested_lambda_values <- value
  validObject(x)
  x
}
#' @rdname lambdas2Test-methods
#' @aliases lambdas2Test
setReplaceMethod("lambdas2Test", signature(x="DNEAobj"),
                 lambdas2TestReplace.DNEAobj)

BICscores.DNEAobj <- function(x){

  x@hyperparameter$BIC_scores
}
#' Access the BIC scores for each lambda value evaluated
#'
#' The function takes as input a \code{\link{DNEAobj}} object and returns
#' the BIC values for each lambda tested during hyperparameter optimization
#' performed via BICtune().
#'
#'
#' @param x A \code{\link{DNEAobj}} object
#' @param value a list containing lists that consist of the liklihood and
#' BIC score for tested lambda values
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{BICtune}}
#' @returns The optimized lambda hyperparameter.
#' @examples
#' #import example data
#' data(dnw)
#'
#' BICscores(dnw)
#' @rdname BICscores-methods
#' @aliases BICscores
#' @export
setMethod("BICscores", signature(x="DNEAobj"),
          BICscores.DNEAobj)

BICscoresReplace.DNEAobj <- function(x, value){

  x@hyperparameter$BIC_scores <- value
  validObject(x)
  x
}
#' @rdname BICscores-methods
#' @aliases BICscores
setReplaceMethod("BICscores", signature(x="DNEAobj"),
                 BICscoresReplace.DNEAobj)

selectionResults.DNEAobj <- function(x){

  x@stable_networks$selection_results
}
#' Access and set the edge selection results from stabilitySelection()
#'
#' The function takes as input a \code{\link{DNEAobj}} object and returns
#' an m x m matrix of selection results for every possible network edge
#' calculated via \code{\link{stabilitySelection}}.
#'
#'
#' @param x A \code{\link{DNEAobj}} object
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{stabilitySelection}},\code{\link{selectionProbabilities}}
#' @returns A \code{\link{DNEAobj}} object after filling the
#'         selection_results section of the stable_networks slot.
#' @examples
#' #import example data
#' data(dnw)
#'
#' selectionResults(dnw)
#' @rdname selectionResults-methods
#' @aliases selectionResults
#' @export
setMethod("selectionResults", signature(x="DNEAobj"),
          selectionResults.DNEAobj)

selectionResultsReplace.DNEAobj <- function(x, value){

  x@stable_networks$selection_results <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("selectionResults", signature(x="DNEAobj"),
                 selectionResultsReplace.DNEAobj)

selectionProbabilities.DNEAobj <- function(x){

  x@stable_networks$selection_probabilities
}
#' Access and set the edge selection probabilities from stabilitySelection()
#'
#' The function takes as input a \code{\link{DNEAobj}} object and returns an
#' m x m matrix of selection probabilities for every possible network
#' edge calculated via \code{\link{stabilitySelection}}.
#'
#'
#' @param x A \code{\link{DNEAobj}} object
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{stabilitySelection}},\code{\link{selectionResults}}
#' @returns A \code{\link{DNEAobj}} object after filling the
#' selection_probabilities section of the stable_networks slot.
#' @examples
#' #import example data
#' data(dnw)
#'
#' selectionProbabilities(dnw)
#' @rdname selectionProbabilities-methods
#' @aliases selectionProbabilities
#' @export
setMethod("selectionProbabilities", signature(x="DNEAobj"),
          selectionProbabilities.DNEAobj)

selectionProbabilitiesReplace.DNEAobj <- function(x, value){

  x@stable_networks$selection_probabilities <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("selectionProbabilities", signature(x="DNEAobj"),
                 selectionProbabilitiesReplace.DNEAobj)

edgeList.DNEAobj <- function(x){

  x@edge_list
}
#' Access the edge list
#'
#' The function takes as input a \code{\link{DNEAobj}} object and returns the
#' edge list created from \code{\link{getNetworks}}.
#'
#'
#' @param x a \code{\link{DNEAobj}} object
#' @param value a data frame of edges in the network
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{getNetworks}},\code{\link{filterNetworks}},
#' \code{\link{getNetworkFiles}}
#' @returns A data frame corresponding to the edge list
#' determined by DNEA
#' @examples
#' #import example data
#' data(dnw)
#'
#' edgeList(dnw)
#' @rdname edgeList-methods
#' @aliases edgeList
#' @export
setMethod("edgeList", signature(x="DNEAobj"),
          edgeList.DNEAobj)

edgeListReplace.DNEAobj <- function(x, value){

  x@edge_list <- value
  validObject(x)
  x
}
#' @rdname edgeList-methods
#' @aliases edgeList
setReplaceMethod("edgeList", signature(x="DNEAobj"),
                 edgeListReplace.DNEAobj)

nodeList.DNEAobj <- function(x){

  x@node_list
}
#' Access the node list
#'
#' The function takes as input a \code{\link{DNEAobj}} object and returns
#' the node list created from \code{\link{createDNEAobject}}.
#'
#'
#' @param x a \code{\link{DNEAobj}} object
#' @param value a data frame of nodes in the network
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}},\code{\link{clusterNet}},
#' \code{\link{getNetworkFiles}}
#' @returns A data frame corresponding to the node
#' list determined by DNEA.
#' @examples
#' #import example data
#' data(dnw)
#'
#' nodeList(dnw)
#' @rdname nodeList-methods
#' @aliases nodeList
#' @export
setMethod("nodeList", signature(x="DNEAobj"),
          nodeList.DNEAobj)

nodeListReplace.DNEAobj <- function(x, value){

  x@node_list <- value
  validObject(x)
  x
}
#' @rdname nodeList-methods
#' @aliases nodeList
setReplaceMethod("nodeList", signature(x="DNEAobj"),
                 nodeListReplace.DNEAobj)

diagnostics.DNEAobj <- function(x){

  x@dataset_summary@diagnostic_values
}
#' Retrieve the diagnostic values for the input expression data
#'
#' This function retrieves teh diagnostic values calculated for
#' the input expression data to \code{\link{createDNEAobject}}
#'
#'
#' @param x \code{\link{DNEAobj}} object or \code{DNEAinputSummary} object
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}}, \code{\link{aggregateFeatures}}
#' @returns Returns the diagnostic values for the input expression data.
#' @examples
#' #import example data
#' data(dnw)
#'
#' diagnostics(dnw)
#' @rdname diagnostics-methods
#' @aliases diagnostics
#' @export
setMethod("diagnostics", signature(x="DNEAobj"),
          diagnostics.DNEAobj)

diagnosticsReplace.DNEAobj <- function(x, value){

  x@dataset_summary$diagnostic_values <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("diagnostics", signature(x="DNEAobj"),
                 diagnosticsReplace.DNEAobj)

diagnostics.DNEAinputSummary <- function(x){

  x@diagnostic_values
}
#' @rdname diagnostics-methods
#' @aliases diagnostics
setMethod("diagnostics", signature(x="DNEAinputSummary"),
          diagnostics.DNEAinputSummary)

diagnosticsReplace.DNEAinputSummary <- function(x, value){

  x@diagnostic_values <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("diagnostics", signature(x="DNEAinputSummary"),
                 diagnosticsReplace.DNEAinputSummary)

datasetSummary.DNEAobj <- function(x){

  x@dataset_summary
}
#' Access the dataset_summary slot of a DNEAobj object
#'
#' This function prints to console the number of samples, number of
#' features, and diagnostic values of the input data to
#' \code{\link{createDNEAobject}} at initiation of the DNEA workflow.
#'
#'
#' @param x A \code{\link{DNEAobj}} object
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}},\code{\link{aggregateFeatures}}
#' @returns The numbers of samples/features and diagnostic values of the
#' input data calculated by \code{\link{createDNEAobject}}.
#' @examples
#' #import example data
#' data(dnw)
#'
#' datasetSummary(dnw)
#' @rdname datasetSummary-methods
#' @aliases datasetSummary
#' @export
setMethod("datasetSummary", signature(x="DNEAobj"),
          datasetSummary.DNEAobj)

datasetSummaryReplace.DNEAobj <- function(x, value){

  x@dataset_summary <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("datasetSummary", signature(x="DNEAobj"),
                 datasetSummaryReplace.DNEAobj)

adjacencyMatrix.DNEAobj <- function(x, weighted=FALSE){

  if(weighted){
    x@adjacency_matrix$weighted_adjacency
  } else if(!weighted){
    x@adjacency_matrix$unweighted_adjacency
  }
}
#' Retrieve the weighted or unweighted adjacency matrix
#'
#' The function takes as input a \code{\link{DNEAobj}} object and returns the
#' weighted or unweighted adjacency matrix determined via
#' \code{\link{getNetworks}}.
#'
#'
#' @param x A \code{\link{DNEAobj}} object
#' @param weighted A boolean indicating whether or not
#' to select the weighted or thresholded (unweighted)
#' adjacency matrix
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{getNetworks}}
#' @returns A matrix corresponding to the adjacency matrix specified.
#' @examples
#' #import example data
#' data(dnw)
#'
#' adjacencyMatrix(dnw, weighted=TRUE)
#' @rdname adjacencyMatrix-methods
#' @aliases adjacencyMatrix
#' @export
setMethod("adjacencyMatrix", signature(x="DNEAobj"),
          adjacencyMatrix.DNEAobj)

adjacencyMatrixReplace.DNEAobj <- function(x, weighted=FALSE, value){

  if(weighted){
    x@adjacency_matrix$weighted_adjacency <- value
  } else if(!weighted){
    x@adjacency_matrix$unweighted_adjacency <- value
  }
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("adjacencyMatrix", signature(x="DNEAobj"),
                 adjacencyMatrixReplace.DNEAobj)

adjacencyGraph.DNEAobj <- function(x, graph){

  x@consensus_clustering@adjacency_graphs[[graph]]
}
#' Retrieve the adjacency graph for the case, control, or joint network
#'
#' The function  returns the adjacency graph made for the case,
#' control, or joint network determined via \code{\link{clusterNet}}.
#'
#'
#' @param x A \code{\link{DNEAobj}} or \code{consensusClusteringResults}
#' object
#' @param graph A character string indicating which of the adjacency
#' graphs to return. Values can be "joint_graph" for the whole graph object,
#' or one of the group values returned by \code{\link{networkGroups}}
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{clusterNet}}
#' @returns An \code{\link{igraph}} graph object corresponding to
#' the specified adjacency graph.
#' @examples
#' #import example data
#' data(dnw)
#'
#' adjacencyGraph(dnw, graph="DM:case")

#' @rdname adjacencyGraph-methods
#' @aliases adjacencyGraph
#' @export
setMethod("adjacencyGraph", signature(x="DNEAobj"),
          adjacencyGraph.DNEAobj)

adjacencyGraphReplace.DNEAobj <- function(x, graph, value){

  x@consensus_clustering@adjacency_graphs$graph <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("adjacencyGraph", signature(x="DNEAobj"),
                 adjacencyGraphReplace.DNEAobj)

adjacencyGraph.consensusClusteringResults <- function(x, graph){

  x@adjacency_graphs$graph
}
#' @rdname adjacencyGraph-methods
#' @aliases adjacencyGraph
setMethod("adjacencyGraph", signature(x="consensusClusteringResults"),
          adjacencyGraph.consensusClusteringResults)

adjacencyGraphReplace.consensusClusteringResults <- function(x, graph, value){

  x@adjacency_graphs$graph <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("adjacencyGraph", signature(x="consensusClusteringResults"),
                 adjacencyGraphReplace.consensusClusteringResults)

summary.consensusClusteringResults <- function(object){

  object@summary
}
#' Retrieve a summary of a consensusClusteringResults object
#'
#' The function takes as input a consensusClusteringResults object and
#' returns a summary of the results of consensus clustering determined
#' via \code{\link{clusterNet}}.
#'
#' @param object A consensusClusteringResults
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{clusterNet}}
#' @returns A data frame that corresponds to a summary of the results
#' of consensus clustering.
#' @rdname summary.consensuClusteringResults-methods
#' @aliases summary.consensusClusteringResults
#' @export
setMethod("summary", signature(object="consensusClusteringResults"),
          summary.consensusClusteringResults)

consensus_clustering.DNEAobj <- function(x){

  x@consensus_clustering
}
#' @keywords internal
#' @noRd
setMethod("consensus_clustering", signature(x="DNEAobj"),
          consensus_clustering.DNEAobj)

CCsummary.DNEAobj <- function(x){

  summary(x@consensus_clustering)
}
#' Retrieves the summary results of consensus clustering
#'
#' The function takes as input a \code{\link{DNEAobj}} object and
#' returns a summary  of the results of consensus clustering determined
#' via \code{\link{clusterNet}}.
#'
#' @param x A \code{\link{DNEAobj}} object
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{clusterNet}}
#' @returns A data frame summary of the consensus clustering
#' results from DNEA.
#' @examples
#' #import example data
#' data(dnw)
#'
#' CCsummary(dnw)
#' @rdname CCsummary-methods
#' @aliases CCsummary
#' @export
setMethod("CCsummary", signature(x="DNEAobj"),
          CCsummary.DNEAobj)

CCsummaryReplace.DNEAobj <- function(x, value){

  x@consensus_clustering@summary <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("CCsummary", signature(x="DNEAobj"),
                 CCsummaryReplace.DNEAobj)

subnetworkMembership.DNEAobj <- function(x){

  x@consensus_clustering@subnetwork_membership
}
#' Retrieve the subnetwork membership for each feature
#'
#' The function takes as input a \code{\link{DNEAobj}} or
#' \code{consensusClusteringResults} object and returns the results of
#' consensus clustering determined via \code{\link{clusterNet}}.
#'
#' @param x A \code{\link{DNEAobj}} object
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{clusterNet}}
#' @returns A data frame that corresponds to the results of
#' consensus clustering.
#' @examples
#' #import example data
#' data(dnw)
#'
#' subnetworkMembership(dnw)
#' @rdname subnetworkMembership-methods
#' @aliases subnetworkMembership
#' @export
setMethod("subnetworkMembership", signature(x="DNEAobj"),
          subnetworkMembership.DNEAobj)

subnetworkMembership.consensusClusteringResults <- function(x){

  x@subnetwork_membership
}
#' @rdname subnetworkMembership-methods
#' @aliases subnetworkMembership
setMethod("subnetworkMembership", signature(x="consensusClusteringResults"),
          subnetworkMembership.consensusClusteringResults)

subnetworkMembershipReplace.DNEAobj <- function(x, value){

  x@consensus_clustering@subnetwork_membership <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("subnetworkMembership", signature(x="DNEAobj"),
                 subnetworkMembershipReplace.DNEAobj)

netGSAresults.DNEAobj <- function(x){

  x@netGSA
}
#' Access the netGSA slot of a DNEAobj object
#'
#' The function takes as input a \code{\link{DNEAobj}} object and returns
#' the netGSA results in the netGSA slot.
#'
#' @param x A \code{\link{DNEAobj}} object
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{runNetGSA}},
#' \code{\link[netgsa:NetGSA]{netgsa::NetGSA()}}
#' @returns A data frame of the results from netGSA.
#' @examples
#' #import example data
#' data(dnw)
#'
#' netGSAresults(dnw)
#' @rdname netGSAresults-methods
#' @aliases netGSAresults
#' @export
setMethod("netGSAresults", signature(x="DNEAobj"),
          netGSAresults.DNEAobj)

netGSAresultsReplace.DNEAobj <- function(x, value){

  x@netGSA <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("netGSAresults", signature(x="DNEAobj"),
                 netGSAresultsReplace.DNEAobj)
