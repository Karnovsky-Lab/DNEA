#' consensusClusteringResults
#'
#' An s4 class to represent the results from consensus
#' clustering within DNEA
#'
#' @slot summary a data frame containing the sub networks
#' as rows and summary information as columns. The columns
#' include: number_of_nodes, number_of_edges,
#' number_of_DE_nodes, and number_of_DE_edges.
#'
#' @slot subnetwork_membership A data frame with the same
#' number of rows as features in the data, and a column
#' indicating which sub network a given feature belongs
#' to, if any.
#'
#' @slot adjacency_graph The resulting adjacency graph
#' from \code{\link{igraph}} created after
#' consensus clustering.
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
#' on the input data to a \code{\link[=DNEA-class]{DNEA}}.
#'
#' @slot num_samples a single-value numeric vector corresponding
#' to the number of samples in the data set.
#'
#' @slot num_features a single-value numeric vector corresponding
#' to the number of features in the data set
#' @slot diagnostic_values a 3x3 data frame with the diagnostic
#' values calculated via \code{\link{createDNEAobject}}.
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

#' DNEA object
#'
#' An s4 class to represent the DNEA workflow
#'
#' @slot project_name A character string name for the experiment.
#'
#' @slot assays A list of matrices, "input_data" being the original
#' non-normalized, non-transformed data, "log_input_data" is the
#' input data log transformed, and "log-scaled_input" is the input
#' data log-transformed and auto-scaled. The row names between the
#' input assays must be identical (the expression data can be
#' accessed via the \code{\link{expressionData}} function).
#' Any other assay input into the DNEA object can be accessed
#' by supplying its name to the assay parameter.
#'
#' @slot metadata A list of information about the data, including a
#' data frame for sample metadata (the row names must match the
#' sample order of the stored expression data), a data frame for
#' feature metadata  (the row names must match the feature order of
#' the stored expression data), a two-level factor corresponding to
#' the two groups in the data, and a character vector the same length as the
#' number of samples corresponding to the group membership for each sample
#' (the user may add additional metadata via the
#' \code{\link{includeMetadata}} function).
#'
#' @slot dataset_summary A \code{DNEAinputSummary} object (can view data
#' via \code{\link{datasetSummary}} and \code{\link{diagnostics}})
#'
#' @slot node_list A data frame containing all of the features in the
#' data set as rows as well as the differential expression analysis
#' results (can view the node list via \code{\link{nodeList}}).
#'
#' @slot edge_list A data frame containing the network edges identified
#' via \code{\link{getNetworks}} (can view the edge list via
#' \code{\link{edgeList}}).
#'
#' @slot hyperparameter A list of results obtained from
#' \code{\link{BICtune}} containing a numeric vector of the
#' lambda values tested during optimization, the resulting
#' Bayesian-information criterion and likelihood scores for each
#' lambda value, and the optimized lambda for analysis (the optimized
#' lambda can be accessed or changed via the
#' \code{\link{optimizedLambda}} function).
#'
#' @slot adjacency_matrix A list of adjacency matrices, one for each
#' experimental condition, jointly estimated via
#' \code{\link{getNetworks}}.
#'
#' @slot stable_networks A list of the selection results and
#' selection probabilities, one for each experimental condition,
#' for every possible feature-feature edge.
#'
#' @slot consensus_clustering A \code{consensusClusteringResults}
#' object containing the results from consensus clustering
#' obtained via the \code{\link{clusterNet}} function.
#'
#' @slot netGSA a data frame containing the results from
#' enrichment analysis performed via \code{\link{runNetGSA}} and
#' the \code{\link[netgsa:NetGSA]{NetGSA}} algorithm. Each row is the
#' results for a given sub network tested for enrichment.
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
#' @returns A DNEA object
#' @inherit createDNEAobject return examples
#'
#' @import methods
#' @name DNEA-class
#' @rdname DNEA-class
#' @aliases DNEA-class
setClass(Class="DNEA",
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

#' collapsed_DNEA
#'
#' An s4 class to represent the DNEA workflow, including collapsing
#' features. This class inherits from the
#' \code{\link[=DNEA-class]{DNEA}} class.
#'
#' @slot original_experiment The \code{\link[=DNEA-class]{DNEA}} object input
#' to \code{\link{aggregateFeatures}}.
#'
#' @slot feature_membership A data frame containing all of the
#' features from the original input data and their corresponding
#' group membership in the new aggregated data.
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{aggregateFeatures}},
#' \code{\link{createDNEAobject}}
#'
#' @returns A collapsed_DNEA object
#' @inherit aggregateFeatures return examples
#'
#' @docType class
#' @import methods
#' @rdname collapsed_DNEA-class
#' @aliases collapsed_DNEA
setClass(Class="collapsed_DNEA",
         slots=c(original_experiment="DNEA",
                 feature_membership="data.frame"),
         contains="DNEA")

#' check validity of "consensusClusteringResults"
#' @docType methods
#' @import methods
#' @noRd
setValidity("consensusClusteringResults", function(object){

  if(ncol(summary(object) != 5)){
    "there was a problem with consensus cluster results object"
  }
  for(i in length(adjacencyGraph(object))){
    if(inherits(adjacencyGraph(object)[i], what="igraph")){
      "There was a problem with adjacency graphs in consensus clustering"
    }
  }
  if(nrow(subnetworkMembership(object)) != nrow(summary(object))){
    "there was a problem with consensus cluster results object"
  }
  if(!all(unique(unlist(subnetworkMembership(object))) %in% c(0, 1))){
    "There was an error in determining sub networks"
  }

  if(sum(summary(object)$number_of_nodes) !=
     ncol(subnetworkMembership(object))){
    "Not all features accounted for in sub networks"
  }
})

#' check validity of "DNEAinputSummary"
#' @docType methods
#' @import methods
#' @noRd
setValidity("DNEAinputSummary", function(object){

  if(!is.numeric(numSamples(object)) | length(numSamples(object)) != 1){
    "there was a problem with dataset summary"
  }
  if(!is.numeric(numFeatures(object)) | length(numFeatures(object)) != 1){
    "there was a problem with dataset summary"
  }
  if(all(dim(diagnostics(object)) != c(3, 2))){
    "There was a problem with diagnostics"
  }
})

#' Check assays slot
#' @aliases DNEA-validator
#' @keywords internal
#' @noRd
assaysCheck <- function(object){
  assays2check <- names(assays.DNEA(object))
  for(i in assays2check){
    if(is.null(assays(object)[[i]])){
      break
    }else if(is.list(assays(object)[[i]])){
      data2check <- assays(object)[[i]]
    }else{data2check <- list(assays(object)[[i]])}
    for(y in seq(length(data2check))){
      if(!(is.matrix(data2check[[y]]))){
        "@assays must be an expression matrix"
      }
      if(length(rownames(data2check[[y]])) !=
         length(unique(rownames(data2check[[y]])))){
        "@assays must be an expression matrix where
        each row is a unique feature."
      }
      if(!(is.numeric(data2check[[y]]))){
        "@assays must be a matrix with numeric values."
      }
      if(numFeatures(object) != nrow(data2check[[y]])){
        "There was a problem with the feature number"
      }
      if(is.null(names(data2check)[y])){
        if(numSamples(object) !=
           ncol(expressionData(x=object, assay=names(assays(object))[1]))){
          "There was a problem with the sample number"
        }

        if(!vector_compare(colnames(data2check[[y]]), sampleNames(object))){
          paste0("The ", i, " matrix should contain the same samples ",
                 "and in the same order as the input data")
        }
      }else if(names(data2check)[y] %in% networkGroups(object)){
        group_samples <- metaData(object, type = "samples")
        samps <- group_samples$conditions == names(data2check)[y]
        group_samples <- rownames(group_samples)[samps]
        if(!vector_compare(colnames(data2check[[y]]), group_samples)){
          paste0("The ", i,": ",  names(data2check)[y],
                 " does not contain the correct group samples")
        }
      }else{
        paste0("The ", i,": ",  names(data2check)[y],
               "does not match a specified experimental group")
      }
      if(all(rownames(data2check[[y]]) !=
             featureNames(object, original=FALSE))){
        "Features are out of order"
      }}}
}
#' Check metadata slot
#' @aliases DNEA-validator
#' @keywords internal
#' @noRd
metadataCheck <- function(object){

  samps <- metaData(object, type="samples")
  if(!(is.data.frame(samps)))"sample metadata should be a data.frame"
  if(!(is.character(samps$samples)))"sample labels should be character strings"

  if(!vector_compare(samps$samples,
                     colnames(expressionData(x=object, assay=names(assays(object))[1])))){
    "sample metadata does not match order of expression data"
  }

  metabs <- metaData(object, type = "features")
  if(!(is.data.frame(metabs)))"feature metadata should be a data.frame"
  if(!(is.character(metabs$feature_names)))"feature names should be characters"
  if(!(is.character(metabs$clean_feature_names))){
    "@metadata$features$clean_Feature_Names should be of class character"
  }
  if(!vector_compare(metabs$clean_feature_names,
                     rownames(expressionData(x=object, assay=names(assays(object))[1])))){
    "feature metadata does not match order of expression data"
  }

  ##check experimental groups data
  if(length(networkGroups(object)) != 2 |
     !is.factor(networkGroupIDs(object))){
    "sample conditions should be a two-level factor"
  }

  if(length(networkGroupIDs(object)) ==
     length(sampleNames(object))){
    if(!vector_compare(names(networkGroupIDs(object)),
                       sampleNames(object))){
      "Group ID's are not aligned with their respective samples"
    }
  }else{
    "There should be only one group label for each sample"
  }
}

#' Check Validity of "DNEA" object
#' @aliases DNEA-validator
#' @docType methods
#' @import methods
#' @noRd
setValidity("DNEA", function(object){

  assaysCheck(object)
  metadataCheck(object)
  validObject(datasetSummary(object))
  validObject(consensus_clustering(object))
  ##check project name
  if(!(is.character(projectName(object)))){
    "@project_name must be a character string"
  }
  ##check nodelist
  metabs <- nodeList(object)$clean_feature_name
  feats <- featureNames(object, original=FALSE)

  if(!vector_compare(metabs[order(metabs)], feats[order(feats)])){
    "Node list features do not match expression data"
  }
  ##check hyperparameter slot
  if(!is.null(BICscores(object)) & !is.null(lambdas2Test(object))){
    if(length(BICscores(object)) != length(lambdas2Test(object))){
      "There was a problem with the tested lambda values"
    }
    Bscores <- unlist(lapply(BICscores(object),
                             function(x) x$BIC))
    if(optimizedLambda(object) !=
       lambdas2Test(object)[match(min(Bscores), Bscores)]){
      "There was a problem with the optimized lambda"
    }}
  ##check adjacency matrices
  if(!is.null(adjacencyMatrix(object, weighted=TRUE))){
    if(all(rownames(adjacencyMatrix(object, weighted=TRUE)) !=
           colnames(expressionData(x=object, assay=names(assays(object))[1])))){
      "there was a problem with the adjacency matrices"
    }}
  ##check stable networks
  if(!is.null(selectionProbabilities(object))){
    for(i in length(selectionProbabilities(object))){
      if(all(dim(selectionProbabilities(object)[[i]]) !=
             c(numFeatures(object), numFeatures(object)))){
        "There was a problem with the feature order for selection probabilites"
      }
      if(any(selectionProbabilities(object)[[i]] > 1) &
         any(selectionProbabilities(object)[[i]] < 0)){
        "There was a problem calculating selection probabilites"
      }}}
})
