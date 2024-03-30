#' consensusClusteringResults
#'
#' An s4 class to represent the results from consensus
#' clustering within DNEA
#'
#' @slot summary a data frame containing the subnetworks as rows and summary
#' information as columns. The columns include: number_of_nodes,
#' number_of_edges, number_of_DE_nodes, number_of_DE_edges
#' @slot subnetwork_membership A data frame with the same number of rows
#' as features in the data, and a column indicating which subnework a given
#' feature belongs to, if any.
#' @slot adjacency_graph The resulting adjacency graph from igraph
#' created after consensus clustering
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{clusterNet}}
#'
#' @returns A consensusClusteringResults object
#' @docType class
#' @import methods
#' @rdname consensusClusteringResults-class
#' @aliases consensusClusteringResults
setClass(Class="consensusClusteringResults",
         slots=c(summary="data.frame",
                 subnetwork_membership="data.frame",
                 adjacency_graphs="list"))

#' DNEAinputSummary
#'
#' An s4 class to represent the results from diagnostic testing
#' on the input data
#'
#' @slot num_samples a single-value numeric vector corresponding to the
#' number of samples in the dataset
#' @slot num_features a single-value numeric vector corresponding to the
#' number of features in the dataset
#' @slot diagnostic_values a 3x3 data frame with the diagnostic values
#' calculated via \code{\link{createDNEAobject}}
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{createDNEAobject}}
#' @returns A DNEAinputSummary object
#' @docType class
#' @import methods
#' @rdname DNEAinputSummary-class
#' @aliases DNEAinputSummary
setClass(Class="DNEAinputSummary",
         slots=c(num_samples="numeric",
                 num_features="numeric",
                 diagnostic_values="data.frame"))

#' DNEAobj
#'
#' An s4 class to represent the DNEA workflow
#'
#' @slot project_name A character string name for the experiment
#'
#' @slot assays A list of matrices, "expression_data" being non-normalized,
#' non-transformed data and "scaled_expression_data" being normalized,
#' log-transformed data. The row names and column names between the two must be
#' identical (can access expression data via \code{\link{expressionData}})
#'
#' @slot metadata A list of information about the data, including a
#' data.frame for sample metadata (rownames must match the sample order of
#' expression data), a data.frame for feature metadata (rownames must match
#' the feature order of expression data), a two-level factor corresponding to
#' the two groups in the data, and a character vector the same length as the
#' number of samples corresponding to the group membership for each sample
#' (the user may add additional metadata via \code{\link{includeMetadata}})
#'
#' @slot dataset_summary A \code{DNEAinputSummary} object (can view data
#' via \code{\link{datasetSummary}} and \code{\link{diagnostics}})
#'
#' @slot node_list A data.frame containing all of the features in the
#' data set as rows as well as the differential expression analysis
#' results (can view the node list via \code{\link{nodeList}})
#'
#' @slot edge_list A data.frame containing the network edges identified
#' via \code{\link{getNetworks}} (can view the edge list via
#' \code{\link{edgeList}})
#'
#' @slot hyperparameter A list of results obtained from
#' \code{\link{BICtune}} containing a numeric vector of the lambda values
#' tested during optimization, the resulting bayesian-information
#' criterion and likelihood scores for each lambda value, and the
#' optimized lambda for analysis (optimized lambda can be accessed or
#' changed via \code{\link{optimizedLambda}})
#'
#' @slot adjacency_matrix A list of adjacency matrices, one for each
#' experimental condition, jointly estimated  via
#' \code{\link{getNetworks}}
#'
#' @slot stable_networks A list of the selection results and
#' selection probabilities, one for each experimental condition,
#' for every possible feature-feature edge
#'
#' @slot consensus_clustering A \code{consensusClusteringResults} object
#' containing the results from consensus clustering via
#' \code{\link{clusterNet}}
#'
#' @slot netGSA a data.frame containing the results from netgsa analysis
#' via \code{\link{runNetGSA}}. Each row is a sub network tested
#' for enrichment
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{expressionData}},\code{\link{includeMetadata}},
#' \code{\link{nodeList}},\code{\link{edgeList}},
#' \code{\link{datasetSummary}},\code{\link{diagnostics}},
#' \code{\link{BICtune}},\code{\link{getNetworks}},
#' \code{\link{stabilitySelection}},\code{\link{clusterNet}},
#' \code{\link{runNetGSA}},\code{\link{selectionProbabilities}},
#' \code{\link{selectionResults}},\code{\link{createDNEAobject}},
#'
#' @returns A DNEAobj object
#' @import methods
#' @name DNEAobj-class
#' @rdname DNEAobj-class
#' @aliases DNEAobj
setClass(Class="DNEAobj",
         slots=c(
           project_name='character',
           assays='list',
           metadata='list',
           dataset_summary='DNEAinputSummary',
           node_list='data.frame',
           edge_list='data.frame',
           hyperparameter='list',
           adjacency_matrix='list',
           stable_networks='list',
           consensus_clustering="consensusClusteringResults",
           netGSA='data.frame')
)

#' collapsed_DNEAobj
#'
#' An s4 class to represent the DNEA workflow, including collapsing
#' features. This class inherits everything from \code{\link{DNEAobj}}
#'
#' @slot original_experiment The DNEAobj object input to
#' \code{\link{aggregateFeatures}}
#' @slot feature_membership A data.frame containing all of the features from
#' the original input data and their corresponding group membership in the
#' collapsed data.
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{aggregateFeatures}},
#' \code{\link{createDNEAobject}}
#'
#' @returns A collapsed_DNEAobj object
#' @docType class
#' @import methods
#' @rdname collapsed_DNEAobj-class
#' @aliases collapsed_DNEAobj
setClass(Class="collapsed_DNEAobj",
         slots=c(original_experiment="DNEAobj",
                 feature_membership="data.frame"),
         contains="DNEAobj")

#' check validity of "consensusClusteringResults"
#' @docType methods
#' @import methods
#' @noRd
setValidity("consensusClusteringResults", function(object){

  if(ncol(object@summary != 5)){
    "there was a problem with consensus cluster results object"
  }
  for(i in length(object@adjacency_graphs)){
    if(inherits(object@adjacency_graphs[i], what="igraph")){
      "There was a problem with adjacency graphs in consensus clustering"
    }
  }
  if(nrow(object@subnetwork_membership) != nrow(object@summary)){
    "there was a problem with consensus cluster results object"
  }
  if(any(object@subnetwork_membership != 1 &
         object@subnetwork_membership != 0)){
    "There was an error in determining metabolic modules"
  }
})

#' check validity of "DNEAinputSummary"
#' @docType methods
#' @import methods
#' @noRd
setValidity("DNEAinputSummary", function(object){

  if(!is.numeric(object@num_samples) | length(object@num_samples) != 1){
    "there was a problem with dataset summary"
  }
  if(!is.numeric(object@num_features) | length(object@num_features) != 1){
    "there was a problem with dataset summary"
  }
  if(all(dim(diagnostics(object)) != c(3, 2))){
    "There was a problem with diagnostics"
  }
})

#' Check Validity of "DNEAobj" object
#' @docType methods
#' @import methods
#' @noRd
setValidity("DNEAobj", function(object){

  ##check project name
  if(!(is.character(projectName(object)))){
    "@project_name must be a character string"
  }

  ##check assays
  for(i in c("input_data", "log_input_data", "scaled_expression_data", "DNEA_scaled_data")){

    if(is.matrix(assays(object)[[i]])){

      data2check <- list(assays(object)[[i]])
    }else if(is.list(assays(object)[[i]])){

      data2check <- assays(object)[[i]]
    }

    for(y in seq(length(data2check))){
      if(!(is.matrix(data2check[[y]]))){
        "@assays must be an expression matrix"
      }
      if(length(rownames(data2check[[y]])) !=
         length(unique(rownames(data2check[[y]])))){
        "@assays must be an expression matrix where each row is a unique feature."
      }
      if(!(is.numeric(data2check[[y]]))){
        "@assays must be a matrix with numeric values."
      }
      if(tryCatch(expr={all(colnames(data2check[[y]]) !=
                            sampleNames(object)[networkGroupIDs(object) == names(data2check)[y]])},
                  warning=function(w){conditionMessage(w)
                    FALSE},
                  error=function(e){conditionMessage(e)
                    FALSE})){
        "Samples are out of order"
      }
      if(all(rownames(data2check[[y]]) !=
             featureNames(object, original=FALSE))){
        "Features are out of order"
      }
    }
  }

  ##check metadata
  if(!(is.data.frame(object@metadata$samples))){
    "@metadata$samples should be of class data.frame"
  }
  if(!(is.character(object@metadata$samples$samples))){
    "@metadata$samples$samples should be of class character"
  }
  if(all(object@metadata$samples$samples !=
         colnames(expressionData(x=object, assay="input_data")))){
    "sample metadata does not match order of expression data"
  }
  if(!(is.data.frame(object@metadata$features))){
    " @metadata$features should be of class data.frame"
  }
  if(!(is.character(object@metadata$features$feature_names))){
    "@metadata$features$feature_names should be of class character"
  }
  if(!(is.character(object@metadata$features$clean_feature_names))){
    "@metadata$features$clean_Feature_Names should be of class character"
  }
  if(all(object@metadata$features$clean_feature_names !=
         rownames(expressionData(x=object, assay="input_data")))){
    "feature metadata does not match order of expression data"
  }
  if(length(networkGroups(object)) != 2){
    "sample conditions should be a two-level factor"
  }
  if(!is.factor(networkGroupIDs(object))){
    "sample conditions should be a factor"
  }
  if(length(networkGroupIDs(object)) !=
     length(sampleNames(object))){
    "There should be only one group label for each sample"
  }

  ##check dataset summary
  if(all(dim(diagnostics(object)) != c(0, 0))){

    if(numSamples(object) !=
       ncol(expressionData(x=object, assay="input_data"))){
      "There was a problem with the dataset summary"
    }
    if(numFeatures(object) !=
       nrow(expressionData(x=object, assay="input_data"))){
      "There was a problem with the dataset summary"
    }
    validObject(object@dataset_summary)
  }

  ##check nodelist
  if(all(dim(nodeList(object)) != c(0, 0))){

    if(!all(nodeList(object)$Features ==
            featureNames(object, original=FALSE))){
      "Node list order does not match expression data"
    }
  }

  ##check edgeList
  if(all(dim(edgeList(object)) != c(0, 0))){

    if(ncol(edgeList(object) != 5)){
      "There was a problem with the edge list"
    }
  }

  ##check hyperparameter slot
  if(!is.null(lambdas2Test(object))){
    if(length(object@hyperparameter$BIC_scores) !=
       length(object@hyperparameter$tested_lambda_values)){
      "There was a problem with the tested lambda values"
    }
    Bscores <- unlist(lapply(object@hyperparameter$BIC_scores, function(x) x$BIC))
    if(optimizedLambda(object) !=
       object@hyperparameter$tested_lambda_values[match(min(Bscores), Bscores)]){
      "There was a problem with the optimized lambda"
    }
  }

  ##check adjacency matrices
  if(!is.null(adjacencyMatrix(object, weighted=TRUE))){

    if(all(rownames(adjacencyMatrix(object, weighted=TRUE)) !=
           rownames(adjacencyMatrix(object, weighted=FALSE)))){
      "there was a problem with the adjacency matrices"
    }
    if(all(colnames(adjacencyMatrix(object, weighted=TRUE)) !=
           colnames(adjacencyMatrix(object, weighted=FALSE)))){
      "there was a problem with the adjacency matrices"
    }
    if(all(rownames(adjacencyMatrix(object, weighted=TRUE)) !=
           colnames(expressionData(x=object, assay="input_data")))){
      "there was a problem with the adjacency matrices"
    }
  }

  ##check stable networks
  if(!is.null(selectionProbabilities(object))){

    for(i in length(selectionProbabilities(object))){

      if(all(dim(selectionProbabilities(object)[[i]]) !=
             c(numFeatures(object), numFeatures(object)))){
        "There was a problem with the feature order for selection probabilites"
      }
      if(any(selectionProbabilities(object)[[i]] > 1) &
         any(selectionProbabilities(object)[[i]] < 0)){
        "There was a problem calculating selectionProbabilites"
      }
    }
  }

  ##check consensus_clustering
  if(!is.null(dim(CCsummary(object)))){

    if(sum(CCsummary(object)$number_of_nodes) !=
       nrow(expressionData(x=object, assay="input_data"))){
      "There was a problem with consensus clustering"
    }

    validObject(object@consensus_clustering)
  }
})
