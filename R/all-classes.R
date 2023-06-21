#'
#'
#'
#' set "igraph" class for consensus cluster
#'
#' @import igraph
#' @keywords internal
#' @noRd
setOldClass('igraph')

#' Set consensusClusteringResults class
#'
#' @import methods
#' @rdname consensusClusteringResults
setClass(Class = "consensusClusteringResults",
         slots = c(summary = "data.frame",
                   subnetwork_membership = "data.frame",
                   adjacency_graphs = "list"))

#' Set DNEAinputData class
#'
#' @import methods
#' @rdname DNEAinputSummary
setClass(Class = "DNEAinputSummary",
         slots = c(num_samples = "numeric",
                   num_features = "numeric",
                   diagnostic_values = "data.frame"))
#' set generic "DNEAresults" class
#'
#'  @import methods
#'  @rdname DNEAresults
setClass(Class = "DNEAresults",
         slots = c(
           project_name = 'character',
           assays = 'list',
           metadata = 'list',
           dataset_summary = 'DNEAinputSummary',
           node_list = 'data.frame',
           edge_list = 'data.frame',
           hyperparameter = 'list',
           adjacency_matrix = 'list',
           stable_networks = 'list',
           joint_graph = 'igraph',
           consensus_clustering = "consensusClusteringResults",
           netGSA = 'list')
)

#' #'Set "collapsed_DNEAresults" class
#' #'
#' #' @import methods
#' #' @noRd
#' setClass(Class = "collapsed_DNEAresults",
#'          slots = c(original = "DNEAresults",
#'                    feature_membership = "list"),
#'          contains = "DNEAobject")
#' #'Set generic "pcorNetwork" class
#' #'
#' #' @import methods
#' #' @noRd
#' setClass(Class = "pcorNetwork",
#'          slots = c(
#'            project_name = 'character',
#'            assays = 'list',
#'            metadata = 'list',
#'            dataset_summary = 'list',
#'            node_list = 'data.frame',
#'            edge_list = 'data.frame',
#'            hyperparameter = 'list',
#'            adjacency_matrix = 'list',
#'            stable_networks = 'list',
#'            joint_graph = 'igraph'
#'          )
#' )
#' #'Set "DNEAobject" class
#' #'
#' #' @import methods
#' #' @noRd
#' setClass(Class = "DNEAobject",
#'          contains = "pcorNetwork",
#'          representation(netGSA = "list",
#'                         consensus_clustering = "list"))


#'Check Validity of "DNEAobject"
#'
#' @import methods
#' @noRd
setValidity("DNEAresults", function(object){
  if(!(is.character(object@project_name))){
    "@project_name must be a character string"
  }
  for (i in length(object@assays)){
    if(!(is.matrix(object@assays[[i]]))){
      "@assays must be an expression matrix"
    }
    if(length(colnames(object@assays[[i]])) != length(unique(colnames(object@assays[[i]])))){
      "@assays must be an expression matrix where each column is a unique feature."
    }
    if(!(is.numeric(object@assays[[i]]))){
      "@assays must be a matrix with numeric values."
    }
    if(all(rownames(object@assays[[i]]) != object@metadata$samples)){
      "Samples are out of order"
    }
    if(all(colnames(object@assays[[i]]) != object@metadata$clean_feature_names)){
      "Features are out of order"
    }
  }
  if(!(is.data.frame(object@metadata$samples))){
    "@metadata$samples should be of class data.frame"
  }
  if(!(is.character(object@metadata$samples$samples))){
    "@metadata$samples$samples should be of class character"
  }
  if(!(is.factor(object@metadata$samples$conditions))){
    "@metadata$samples$conditions should be of class factor"
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
})

############################################################################



