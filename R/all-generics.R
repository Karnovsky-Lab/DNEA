#' @include all-classes.R
NULL

#' @rdname BICtune-methods
#' @aliases BICtune
#' @export
setGeneric("BICtune", function(object,
                               lambda_values,
                               interval=1e-3,
                               informed=TRUE,
                               eps_threshold=1e-6,
                               eta_value=0.1,
                               BPPARAM=bpparam(),
                               BPOPTIONS=bpoptions())
  standardGeneric("BICtune"))

#' @rdname projectName-methods
#' @aliases projectName
#' @export
setGeneric("projectName", function(x) standardGeneric("projectName"))

#' @rdname expressionData-methods
#' @aliases expressionData
#' @export
setGeneric("expressionData", function(x, assay)
  standardGeneric("expressionData"))

#' @rdname assays-methods
#' @aliases assays
#' @keywords internal
#' @noRd
setGeneric("assays", function(x)
  standardGeneric("assays"))

#' @rdname networkGroups-methods
#' @aliases networkGroups
#' @export
setGeneric("networkGroups", function(x)
  standardGeneric("networkGroups"))

#' @rdname networkGroupIDs-methods
#' @aliases networkGroupIDs
#' @export
setGeneric("networkGroupIDs", function(x)
  standardGeneric("networkGroupIDs"))

#' @rdname networkGroupIDs-methods
#' @aliases networkGroupIDs
#' @export
setGeneric("networkGroupIDs<-", function(x, value)
  standardGeneric("networkGroupIDs<-"))

#' @rdname metaData-methods
#' @aliases metaData includeMetadata
#' @export
setGeneric("metaData", function(x, type)
  standardGeneric("metaData"))

#' @rdname metaData-methods
#' @aliases metaData includeMetadata
#' @export
setGeneric("metaData<-", function(x, type, value)
  standardGeneric("metaData<-"))

#' @rdname sampleNames-methods
#' @aliases sampleNames
#' @export
setGeneric("sampleNames", function(x)
  standardGeneric("sampleNames"))

#' @rdname featureNames-methods
#' @aliases featureNames
#' @export
setGeneric("featureNames", function(x, original=TRUE)
  standardGeneric("featureNames"))

#' @rdname numFeatures-methods
#' @aliases numFeatures
#' @export
setGeneric("numFeatures", function(x)
  standardGeneric("numFeatures"))

#' @rdname numSamples-methods
#' @aliases numSamples
#' @export
setGeneric("numSamples", function(x)
  standardGeneric("numSamples"))

#' @rdname optimizedLambda-methods
#' @aliases optimizedLambda
#' @export
setGeneric("optimizedLambda", function(x)
  standardGeneric("optimizedLambda"))

#' @rdname optimizedLambda-methods
#' @aliases optimizedLambda
#' @export
setGeneric("optimizedLambda<-", function(x, value)
  standardGeneric("optimizedLambda<-"))

#' @rdname lambdas2Test-methods
#' @aliases lambdas2Test
#' @export
setGeneric("lambdas2Test", function(x)
  standardGeneric("lambdas2Test"))

#' @rdname lambdas2Test-methods
#' @aliases lambdas2Test
#' @export
setGeneric("lambdas2Test<-", function(x, value)
  standardGeneric("lambdas2Test<-"))

#' @rdname BICscores-methods
#' @aliases BICscores
#' @export
setGeneric("BICscores", function(x)
  standardGeneric("BICscores"))

#' @rdname BICscores-methods
#' @aliases BICscores
#' @export
setGeneric("BICscores<-", function(x, value)
  standardGeneric("BICscores<-"))

#' @rdname selectionResults-methods
#' @aliases selectionResults
setGeneric("selectionResults", function(x)
  standardGeneric("selectionResults"))

#' @keywords internal
#' @noRd
setGeneric("selectionResults<-", function(x, value)
  standardGeneric("selectionResults<-"))

#' @rdname selectionProbabilities-methods
#' @aliases selectionProbabilities
setGeneric("selectionProbabilities", function(x)
  standardGeneric("selectionProbabilities"))

#' @keywords internal
#' @noRd
setGeneric("selectionProbabilities<-", function(x, value)
  standardGeneric("selectionProbabilities<-"))

#' @rdname nodeList-methods
#' @aliases nodeList
#' @export
setGeneric("nodeList", function(x)
  standardGeneric("nodeList"))

#' @rdname nodeList-methods
#' @aliases nodeList
setGeneric("nodeList<-", function(x, value)
  standardGeneric("nodeList<-"))

#' @rdname edgeList-methods
#' @aliases edgeList
#' @export
setGeneric("edgeList", function(x)
  standardGeneric("edgeList"))

#' @rdname edgeList-methods
#' @aliases edgeList
setGeneric("edgeList<-", function(x, value)
  standardGeneric("edgeList<-"))

#' @rdname diagnostics-methods
#' @aliases diagnostics
setGeneric("diagnostics", function(x)
  standardGeneric("diagnostics"))

#' @keywords internal
#' @noRd
setGeneric("diagnostics<-", function(x, value)
  standardGeneric("diagnostics<-"))

#' @rdname datasetSummary-methods
#' @aliases datasetSummary
setGeneric("datasetSummary", function(x)
  standardGeneric("datasetSummary"))

#' @keywords internal
#' @noRd
setGeneric("datasetSummary<-", function(x, value)
  standardGeneric("datasetSummary<-"))

#' @rdname adjacencyMatrix-methods
#' @aliases adjacencyMatrix
setGeneric("adjacencyMatrix", function(x, weighted)
  standardGeneric("adjacencyMatrix"))

#' @keywords internal
#' @noRd
setGeneric("adjacencyMatrix<-", function(x, weighted, value)
  standardGeneric("adjacencyMatrix<-"))


#' @rdname subnetworkMembership-methods
#' @aliases subnetworkMembership
#' @export
setGeneric("subnetworkMembership", function(x)
  standardGeneric("subnetworkMembership"))

#' @keywords internal
#' @noRd
setGeneric("subnetworkMembership<-", function(x, value)
  standardGeneric("subnetworkMembership<-"))

#' @rdname adjacencyGraph-methods
#' @aliases adjacencyGraph
#' @export
setGeneric("adjacencyGraph", function(x, graph)
  standardGeneric("adjacencyGraph"))

#' @keywords internal
#' @noRd
setGeneric("adjacencyGraph<-", function(x, graph, value)
  standardGeneric("adjacencyGraph<-"))

#' @rdname CCsummary-methods
#' @aliases CCsummary
#' @export
setGeneric("CCsummary", function(x)
  standardGeneric("CCsummary"))

#' @keywords internal
#' @noRd
setGeneric("CCsummary<-", function(x, value)
  standardGeneric("CCsummary<-"))

#' @rdname netGSAresults-methods
#' @aliases netGSAresults
#' @export
setGeneric("netGSAresults", function(x)
  standardGeneric("netGSAresults"))


#' @keywords internal
#' @noRd
setGeneric("netGSAresults<-", function(x, value)
  standardGeneric("netGSAresults<-"))


#' @rdname filterNetworks-methods
#' @aliases filterNetworks
#' @export
setGeneric("filterNetworks", function(data, pcor, top_percent_edges)
  standardGeneric("filterNetworks"))
