
#'
#' @include all-classes.R
#' @rdname expressionData
#' @export
setGeneric("expressionData", function(x, type) standardGeneric("expressionData"))

#'
#' @rdname networkGroups
#' @export
setGeneric("networkGroups", function(x) standardGeneric("networkGroups"))

#'
#' @rdname networkGroups
#' @export
setGeneric("networkGroups<-", function(x, value) standardGeneric("networkGroups<-"))

#'
#' @rdname networkGroupIDs
#' @export
setGeneric("networkGroupIDs", function(x) standardGeneric("networkGroupIDs"))

#'
#' @rdname sampleNames
#' @export
setGeneric("sampleNames", function(x) standardGeneric("sampleNames"))

#'
#' @rdname featureNames
#' @export
setGeneric("featureNames", function(x, original = TRUE) standardGeneric("featureNames"))

#'
#' @rdname numFeatures
#' @export
setGeneric("numFeatures", function(x) standardGeneric("numFeatures"))

#'
#' @rdname numSamples
#' @export
setGeneric("numSamples", function(x) standardGeneric("numSamples"))

#'
#' @rdname optimizedLambda
#' @export
setGeneric("optimizedLambda", function(x) standardGeneric("optimizedLambda"))

#'
#' @rdname optimizedLambda
#' @export
setGeneric("optimizedLambda<-", function(x, value) standardGeneric("optimizedLambda<-"))

#'
#' @rdname lambdas2Test
#' @export
setGeneric("lambdas2Test", function(x) standardGeneric("lambdas2Test"))

#'
#' @rdname lambdas2Test
#' @export
setGeneric("lambdas2Test<-", function(x, value) standardGeneric("lambdas2Test<-"))

#'
#' @rdname BICscores
#' @export
setGeneric("BICscores", function(x) standardGeneric("BICscores"))

#'
#' @rdname BICscores
#' @export
setGeneric("BICscores<-", function(x, value) standardGeneric("BICscores<-"))

#'
#' @rdname selectionResults
#' @keywords internal
setGeneric("selectionResults", function(x) standardGeneric("selectionResults"))

#'
#' @rdname selectionResults
#' @keywords internal
setGeneric("selectionResults<-", function(x, value) standardGeneric("selectionResults<-"))

#'
#' @rdname selectionProbabilities
#' @keywords internal
setGeneric("selectionProbabilities", function(x) standardGeneric("selectionProbabilities"))

#'
#' @rdname selectionProbabilities
#' @keywords internal
setGeneric("selectionProbabilities<-", function(x, value) standardGeneric("selectionProbabilities<-"))
#'
#' @rdname nodeList
#' @export
setGeneric("nodeList", function(x) standardGeneric("nodeList"))

#'
#' @rdname nodeList
#' @keywords internal
setGeneric("nodeList<-", function(x, value) standardGeneric("nodeList<-"))

#'
#' @rdname edgeList
#' @export
setGeneric("edgeList", function(x) standardGeneric("edgeList"))

#'
#' @rdname edgeList
#' @keywords internal
setGeneric("edgeList<-", function(x, value) standardGeneric("edgeList<-"))

#'
#' @rdname diagnostics
#' @keywords internal
setGeneric("diagnostics", function(x) standardGeneric("diagnostics"))

#'
#' @rdname diagnostics
#' @keywords internal
setGeneric("diagnostics<-", function(x, value) standardGeneric("diagnostics<-"))

#'
#' @rdname datasetSummary
#' @keywords internal
setGeneric("datasetSummary", function(x) standardGeneric("datasetSummary"))

#'
#' @rdname datasetSummary
#' @keywords internal
setGeneric("datasetSummary<-", function(x, value) standardGeneric("datasetSummary<-"))

#'
#' @rdname adjacencyMatrix
#' @keywords internal
setGeneric("adjacencyMatrix", function(x, weighted) standardGeneric("adjacencyMatrix"))

#'
#' @rdname adjacencyMatrix
#' @keywords internal
setGeneric("adjacencyMatrix<-", function(x, weighted, value) standardGeneric("adjacencyMatrix<-"))


#' @rdname subnetworkMembership
#' @export
setGeneric("subnetworkMembership", function(x) standardGeneric("subnetworkMembership"))

#' @rdname adjacencyGraph
#' @export
setGeneric("adjacencyGraph", function(x, graph) standardGeneric("adjacencyGraph"))

#'
#' @rdname adjacencyGraph
#' @export
setGeneric("adjacencyGraph<-", function(x, graph, value) standardGeneric("adjacencyGraph<-"))

#'
#' @rdname CCsummary
#' @export
setGeneric("CCsummary", function(x) standardGeneric("CCsummary"))

#'
#' @rdname netGSAresults
#' @export
setGeneric("netGSAresults", function(x) standardGeneric("netGSAresults"))

#'
#' @rdname netGSAresults
#' @export
setGeneric("netGSAresults<-", function(x, value) standardGeneric("netGSAresults<-"))

#'
#' @rdname filterNetworks
#' @export
setGeneric("filterNetworks", function(data, pcor, top_percent_edges) standardGeneric("filterNetworks"))
