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
      "  Samples -  ", paste0('There are ',length(sampleNames(object)), ' samples.'), "\n",
      "  Features -  ", paste0('There are ',length(featureNames(object)), ' Features.'), "\n",
      "  Conditions -  ", paste0("control: ", networkGroups(object)[[1]],
                                 " case: ", networkGroups(object)[[2]]), "\n",
      "  Optimized Lambda - ", paste0("lambda: ",optimizedLambda(object),
                                      " will be used in analysis."), "\n",
      "  Sub-clusters: ", paste0("there are ", length(unique(nodeList(object)[["membership"]])),
                                 " sub-clusters."),
      sep = ""
  )
})

#'#' Show function will display general information about the data present in the DNEAinputSummary slots
#' @param object A DNEAinputSummary object
#' @export
#' @noRd
setMethod("show", "DNEAinputSummary", function(object) {
  cat(paste0(is(object)[[1]], "\n",
             "  Number of Samples  -  ", numSamples(object), "\n",
             "  Number of Features  -  ", numFeatures(object), "\n"),
      sep = "")
  print(as.matrix(diagnostics(object)))
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

networkGroupIDs.DNEAresults <- function(x){
  x@metadata$network_group_IDs
}
#' networkGroupIDs returns the group information for samples in a DNEAobject,
#'  or DNEAobject_collapsed
#'
#' This function takes in a DNEAobject or DNEAobject_collapsed and return the unique condition
#'  labels located in the metadata slot of the object. You can also change the sample
#'  variable being used to determine groups.
#'
#'
#' @param x A DNEAobject or DNEAobject_collapsed
#' @return A vector of the unique condition labels
#'
#' @rdname networkGroupIDs
#' @export
setMethod("networkGroupIDs", signature(x = "DNEAresults"), networkGroupIDs.DNEAresults)

#' projectName() returns the name of the current experiment
#'
#' projectName() returns the name of the current experiment and allows the user to change the name as well
#'
#' @param x A DNEAresults object
#' @param value A string that corresponds to the new project name
#' @return The name of the DNEA experiment
#'
#' @rdname projectName
#' @export
setMethod("projectName", signature(x= "DNEAresults"), function(x){

  x@project_name
})

#'
#' @rdname projectName
#' @export
setReplaceMethod("projectName", signature(x= "DNEAresults"), function(x, value){

  x@project_name <- value
  validObject(x)
  x
})
#' networkGroups retrieves the condition values for each sample from the metadata slot.
#'
#' This function takes in a DNEAobject, or DNEAobject_collapsed object and
#'  return the condition values located in the metatdata slot of the object
#'
#' @param x A DNEAobject, or DNEAobject_collapsed object
#' @param value A vector of group names that specifies the condition to test the data on. It should have exactly two unique
#' groups and be exactly as long as the number of samples in the data.
#' @return A vector of the condition values
#'
#' @rdname networkGroups
#' @export
setMethod("networkGroups", signature(x = "DNEAresults"), function(x){

  x@metadata$network_groups
})

#' @rdname networkGroups
#' @export
setMethod("networkGroups<-", signature(x = "DNEAresults"), function(x, value){

  #set new group ID's
  if(length(unique(x@metadata$samples[[value]])) >2) stop("DNEA requires exactly two groups!")
  x@metadata$network_group_IDs <- x@metadata$samples[[value]]

  #update group levels
  x@metadata$network_groups <- levels(x@metadata$samples[[value]])
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
  x@dataset_summary@num_features
}
#' numFeatures finds the total number of features in the dataset
#'
#' @param x A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return The number of features in the dataset
#'
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
#' numSamples finds the total number of samples in the dataset
#'
#' @param x A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return The number of samples in the dataset
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
#' optimizedLambda returns the lambda value used in analysis.
#'
#' The function takes as input a pcorNetwork, DNEAobject, or DNEAobject_collapsed object and returns
#'  the hyperparameter (lambda) that is currently being used for the analysis
#'
#' @param x A pcorNetwork, DNEAobjct, or DNEAobject_collapsed object
#' @param value A numeric value between 0 and 1 corresponding to the lambda parameter to use in the analysis
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
#' @param value A vector of lambda values to optimize in BICtune()
#'
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
diagnostics.DNEAresults <- function(x){

  x@dataset_summary@diagnostic_values
}
#' diagnostics retrieves the diagnostic values for the input expression data
#'
#' FOR INTERNAL USE ONLY - The function takes as input a pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' and returns the diagnostic values for the input expression data
#'
#' @param a pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @return returns the diagnostic values for the input expression data
#'
#' @rdname diagnostics
#' @keywords internal
setMethod("diagnostics", signature(x = "DNEAresults"), diagnostics.DNEAresults)

#' @rdname diagnostics
#' @keywords internal
setReplaceMethod("diagnostics", signature(x = "DNEAresults"), function(x, value){

  x@dataset_summary$diagnostic_values <- value
  validObject(x)
  x
})

#' @rdname diagnostics
#' @export
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

#' datasetSummary is a setter/getter function for the dataset_summary slot of a DNEAresults object
#'
#' @param x A DNEA results object
#' @returns The DNEAinputSummary object if retrieving data, or the input object after populating the dataset_summary
#' slot
#'
#' @rdname datasetSummary
#' @keywords internal
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

#' summary retrieves the summary results of consensus clustering
#'
#' The function takes as input a consensusClusteringResults object and returns a summary of the results of
#' consensus clustering determined via runConsensusCluster().
#'
#' @param object A consensusClusteringResults
#' @return A data.frame that corresponds to a summary of the results of consensus clustering
#'
#' @rdname summary
#' @export
setMethod("summary", signature(object = "consensusClusteringResults"), function(object){

  object@summary
})

#' CCsummary retrieves the summary results of consensus clustering
#'
#' The function takes as input a pcorNetwork, DNEAobject, or DNEAobject_collapsed object and returns a summary
#' of the results of consensus clustering determined via runConsensusCluster().
#'
#' @param x A pcorNetwork, DNEAobject, or DNEAobject_collapsed
#' @returns A data.frame summary of the consensus clustering results from DNEA
#'
#' @rdname CCsummary
#' @export
setMethod("CCsummary", signature(x = "DNEAresults"), function(x){

  summary(x@consensus_clustering)
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

#' netGSAresults accesses the netGSA slot of a DNEAresults object
#'
#' The function takes as input a pcorNetwork, DNEAobject, DNEAobject_collapsed, or consensusClusteringResults object
#' and returns the netGSA results.
#' **FOR INTERNAL USE ONLY - The results can also be inserted into the object via the setter functionality of the function.**
#'
#' @param x A DNEAresults object
#' @returns The results from netGSA
#'
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

filterNetworks.DNEAresults <- function(data,
                                       pcor,
                                       top_percent_edges){

  ##grab adjacency matrices
  weighted_adjacency_matrices <- adjacencyMatrix(data, weighted = TRUE)

  ##start new unweighted adjacency list
  unweighted_adjacency_matrices <- vector("list", length = length(weighted_adjacency_matrices))
  names(unweighted_adjacency_matrices) <- names(weighted_adjacency_matrices)

  ##filter networks
  #can't provide both filters
  if(!missing(pcor) & !missing(top_percent_edges)){

    stop("Only pcor or top_percent_edges should be specified, not both!")
  }else if(!missing(pcor)){

    #filter by pcor
    for(k in names(weighted_adjacency_matrices)){
      weighted_adjacency_matrices[[k]][abs(weighted_adjacency_matrices[[k]]) < pcor] <- 0
      unweighted_adjacency_matrices[[k]] <- weighted_adjacency_matrices[[k]] != 0
    }
  }else if(!missing(top_percent_edges)){

    #filter networks to top x% strongest edges
    for(k in names(weighted_adjacency_matrices)){

      #grab quantile value
      percent_cutoff <- quantile(abs(weighted_adjacency_matrices[[k]][weighted_adjacency_matrices[[k]] != 0]),
                                 probs = (1 - top_percent_edges))

      #filter
      weighted_adjacency_matrices[[k]][abs(weighted_adjacency_matrices[[k]]) < percent_cutoff] <- 0
      unweighted_adjacency_matrices[[k]] <- weighted_adjacency_matrices[[k]] != 0
    }
  }else{

    #must provide one filter
    stop("Neither pcor nor top_percent_edges were specified - No filtering was performed!")
  }

  #store the adjacency matrices in DNEAresults object
  adjacencyMatrix(x = data, weighted = TRUE) <- weighted_adjacency_matrices
  adjacencyMatrix(x = data, weighted = FALSE) <- unweighted_adjacency_matrices

  ##update edge list
  #initiate output dataframe
  pairs <- combn(as.character(featureNames(data)), 2, simplify=FALSE)
  edge_list <- data.frame(Metabolite.A=rep(0,length(pairs)), Metabolite.B=rep(0,length(pairs)),
                          pcor.0=rep(0,length(pairs)), pcor.1=rep(0,length(pairs)),
                          check.names = FALSE)

  #concatenate results into dataframe
  edge_list[,1:2] <- do.call(rbind, pairs)
  edge_list[,3] <- lowerTriangle(weighted_adjacency_matrices[[1]])
  edge_list[,4] <- lowerTriangle(weighted_adjacency_matrices[[2]])
  edge_list$edge <- rep(NA, length(pairs)) #non-edge
  edge_list$edge[edge_list$pcor.0 != 0 & edge_list$pcor.1 == 0] <- names(weighted_adjacency_matrices)[[1]] #control
  edge_list$edge[edge_list$pcor.0 == 0 & edge_list$pcor.1 != 0] <- names(weighted_adjacency_matrices)[[2]] #case
  edge_list$edge[edge_list$pcor.0 != 0 & edge_list$pcor.1 != 0] <- "Both" #Both
  edge_list <- edge_list[!is.na(edge_list$edge),] #remove non-edges
  rownames(edge_list) <- NULL

  #replace edgelist
  edgeList(data) <- edge_list


  #print message for total edges
  message(paste0(names(unweighted_adjacency_matrices)[[1]]," network specific edges: ", sum(unweighted_adjacency_matrices[[1]])/2), appendLF = TRUE)
  message(paste0(names(unweighted_adjacency_matrices)[[2]]," network specific edges: ", sum(unweighted_adjacency_matrices[[2]])/2), appendLF = TRUE)
  message(rep("-", 35), appendLF = TRUE)
  message(paste0("Number of edges shared by both networks: ", sum(edge_list$edge == "Both")))
  message(paste0("Total number of edges in dataset: ", nrow(edge_list)))

  return(data)
}
#' filterNetworks() will reduce the adjacency matrices to only the edges that meet the filter conditions
#'
#' filterNetworks() takes as input a DNEAresults object and allows the user to filter the network edges by one of two ways:
#'
#' 1. The networks can be filtered to only include edges greater than or equal to a specified partial correlation (pcor) value.
#'
#' 2. The networks can be filtered to only include the strongest x% of edges by pcor value.
#'
#' Filtering is performed on the case and control adjacency matrix separately.
#'
#' @param data A DNEAresults object
#' @param pcor A pcor value of which to threshold the adjacency matrices. Edges with pcor values <= to this value will be removed.
#' @param top_percent_edges A value between 0-1 that corresponds to the top x% edges to keep in the networks
#' ie. top_percent_edges = 0.1 will keep only the top 10% strongest edges in the networks.
#'
#' @returns The input object after filtering the egdes in the network according to the specified parameters
#'
#' @importFrom stats quantile
#' @rdname filterNetworks
#' @export
setMethod("filterNetworks", signature(data = "DNEAresults"), filterNetworks.DNEAresults)

filterNetworks.list <- function(data, pcor, top_percent_edges){

  ##rename input variable
  weighted_adjacency_matrices <- data

  ##start new unweighted adjacency list
  unweighted_adjacency_matrices <- vector("list", length = length(weighted_adjacency_matrices))
  names(unweighted_adjacency_matrices) <- names(weighted_adjacency_matrices)

  ##filter networks
  #can't provide both filters
  if(!missing(pcor) & !missing(top_percent_edges)){

    stop("Only pcor or top_percent_edges should be specified, not both!")
  }else if(!missing(pcor)){

    #filter by pcor
    for(k in names(weighted_adjacency_matrices)){
      weighted_adjacency_matrices[[k]][abs(weighted_adjacency_matrices[[k]]) < pcor] <- 0
      unweighted_adjacency_matrices[[k]] <- weighted_adjacency_matrices[[k]] != 0
    }
  }else if(!missing(top_percent_edges)){

    #filter networks to top x% strongest edges
    for(k in names(weighted_adjacency_matrices)){

      #grab quantile value
      percent_cutoff <- quantile(abs(weighted_adjacency_matrices[[k]][weighted_adjacency_matrices[[k]] != 0]),
                                 probs = (1 - top_percent_edges))

      #filter
      weighted_adjacency_matrices[[k]][abs(weighted_adjacency_matrices[[k]]) < percent_cutoff] <- 0
      unweighted_adjacency_matrices[[k]] <- weighted_adjacency_matrices[[k]] != 0
    }
  }else{

    #must provide one filter
    stop("Neither pcor nor top_percent_edges were specified - No filtering was performed!")
  }

  ##print message for total edges
  message(paste0("Number of edges in ", names(unweighted_adjacency_matrices)[[1]],": ", sum(unweighted_adjacency_matrices[[1]])/2), appendLF = TRUE)
  message(paste0("Number of edges in ", names(unweighted_adjacency_matrices)[[2]],": ", sum(unweighted_adjacency_matrices[[2]])/2), appendLF = TRUE)

  return(list(weighted = weighted_adjacency_matrices, unweighted = unweighted_adjacency_matrices))
}
#'
#' @rdname filterNetworks
#' @keywords internal
setMethod("filterNetworks", signature(data = "list"), filterNetworks.list)

