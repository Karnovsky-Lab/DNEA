
#' set "igraph" class for consensus cluster
#'
#' @import igraph
#' @keywords internal
#' @noRd
setOldClass('igraph')

#'Set generic "pcorNetwork" class
#'
#' @import methods
#' @noRd
setClass(Class = "pcorNetwork",
         slots = c(
           project_name = 'character',
           assays = 'list',
           metadata = 'list',
           dataset_summary = 'list',
           node_list = 'data.frame',
           edge_list = 'data.frame',
           hyperparameter = 'list',
           adjacency_matrix = 'list',
           stable_networks = 'list',
           joint_graph = 'igraph'
         )
)
#'Set "DNEAobject" class
#'
#' @import methods
#' @noRd
setClass(Class = "DNEAobject",
         representation(netGSA_results = "list"),
         contains = "pcorNetwork")
# setClass(Class = "DNEAobject",
#          slots = c(
#            project_name = 'character',
#            assays = 'list',
#            metadata = 'list',
#            dataset_summary = 'list',
#            node_list = 'data.frame',
#            edge_list = 'data.frame',
#            hyperparameter = 'list',
#            adjacency_matrix = 'list',
#            stable_networks = 'list',
#            joint_graph = 'igraph',
#            netGSA_results = 'list',
#            feature_membership = 'list'
#          )
# )
#'Check Validity of "DNEAobject"
#'
#' @import methods
#' @noRd
setValidity("DNEAobject", function(object){
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

#' Show function will display general information about the data present in the DNEAobject slots
#' @param object A DNEAobject
#' @export
#' @noRd
setMethod("show", "DNEAobject", function(object) {
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

#' split_by_condition splits the input data by condition
#'
#' split_by_condition will separate the input expression matrix into a list of matrices. Each matrix
#' corresponds to the expression data for one condition specified by condition_levels
#'
#' @param dat a matrix of expression data wherein the samples are rows and features are columns.
#' @param condition_levels A list or vector of the unique conditions present in dat
#' @param condition_by_sample A list or vector of the condition value corresponding to each
#'        sample.
#' @return A list of expression matrices
#' @keywords internal
split_by_condition <- function(dat, condition_levels, condition_by_sample){

  #create key for separating the data by key and running diagnostic tests, feature DE calculations
  separated_conditions_data <- vector(mode = 'list', length = length(condition_levels))
  names(separated_conditions_data) <- condition_levels

    for(cond in condition_levels){
      separated_conditions_data[[cond]] <- t(dat)[,(condition_by_sample == cond)]

    }
  return(separated_conditions_data)

}
#' expressionData retrieves un-scaled expression data from the assays slot.
#'
#' This function takes in a DNEAobject and return the un-scaled expression matrix located in the
#' assays slot of the object. You can also use it to input expression data into the DNEAobject
#'
#' @param x A DNEAobject
#' @param type "input" will return the original data, and "normalized" will return the scaled data
#' @return The indicated expression expression matrix
#' @export
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
#'
setGeneric("expressionData", function(x, type) standardGeneric("expressionData"))
setMethod("expressionData",signature(x = "pcorNetwork"), expressionData.pcorNetwork)


# setGeneric("expressionData<-",
#            function(x, value) standardGeneric("expressionData<-"))
# setMethod("expressionData<-", "DNEAobject", function(x, value) {
#   x@assays[['expresion_data']] <- value
#   validObject(x)
#   return(x)
# })

#' conditionLevels retrieves the unique condition names for the dataset
#'
#' This function takes in a DNEAobject and return the unique condition labels located in the
#' dataset_summary slot of the object
#'
#' @param x A DNEAobject
#' @return A vector of the unique condition labels
#' @export
conditionLevels.pcorNetwork <- function(x){
  x@dataset_summary$condition_levels
}

setGeneric("conditionLevels", function(x) standardGeneric("conditionLevels"))
setMethod("conditionLevels", signature(x = "pcorNetwork"), conditionLevels.pcorNetwork)

#' condition retrieves the condition values for each sample from the metadata slot.
#'
#' This function takes in a DNEAobject and return the condition values located in the
#' metatdata slot of the object
#'
#' @param x A DNEAobject
#' @return A vector of the condition values
#' @export
setGeneric("conditions",
           function(x) standardGeneric("conditions"))
setMethod("conditions", signature(x = "pcorNetwork"), function(x) x@metadata$samples$conditions)
setGeneric("conditions<-", function(x, value) standardGeneric("conditions<-"))
setMethod("conditions<-", signature(x = "pcorNetwork"), function(x, value){
  x@metadata$samples$conditions <- x@metadata$samples[[value]]
})

#' sampleNames retrieves the sample names from the metadata slot.
#'
#' This function takes in a DNEAobject and return the sample names located in the
#' metadata slot of the object.
#'
#' @param x A DNEAobject
#' @return A vector of sample names
#' @keywords internal
sampleNames.pcorNetwork <- function(x){
  x@metadata[["samples"]]$samples
}
setGeneric("sampleNames",
           function(x) standardGeneric("sampleNames"))
setMethod("sampleNames", signature(x = "pcorNetwork"), sampleNames.pcorNetwork)

#' featureNames retrieves the feature names from the metadata slot.
#'
#' This function takes in a DNEAobject and return the feature names located in the
#' metadata slot of the object.
#'
#' @param x A DNEAobject
#' @return A vector of feature names
#' @keywords internal
featureNames.pcorNetwork <- function(x){
  x@metadata[["features"]]$feature_names
}
setGeneric("featureNames",
           function(x) standardGeneric("featureNames"))
setMethod("featureNames", signature(x = "pcorNetwork"), featureNames.pcorNetwork)

#' cleanFeatureNames retrieves the feature names from clean_feature names in the metadata slot.
#'
#' This function takes in a DNEAobject and return the feature names modified to avoid errors with
#' native R clashes. The names are located in the metadata slot of the object.
#'
#' @param x A DNEAobject
#' @return A vector of modified feature names
#' @keywords internal
cleanFeatureNames.pcorNetwork <- function(x){
  x@metadata[["features"]]$clean_feature_names
}
setGeneric("cleanFeatureNames", function(x) standardGeneric("cleanFeatureNames"))
setMethod("cleanFeatureNames", signature(x = "pcorNetwork"), cleanFeatureNames.pcorNetwork)

#' numFeatures finds the total number of features in the dataset
#'
#' @param x A DNEAobject
#' @return The number of features in the dataset
#' @keywords internal
numFeatures.pcorNetwork <- function(x){
  x@dataset_summary$num_features
}
setGeneric("numFeatures",
           function(x) standardGeneric("numFeatures"))
setMethod("numFeatures", signature(x = "pcorNetwork"), numFeatures.pcorNetwork)

#' numSamples finds the total number of samples in the dataset
#'
#' @param x A DNEAobject
#' @return The number of samples in the dataset
#' @keywords internal
numSamples.pcorNetwork <- function(x){
  x@dataset_summary$num_samples
}
setGeneric("numSamples",
           function(x) standardGeneric("numSamples"))
setMethod("numSamples", signature(x = "pcorNetwork"), function(x) numSamples.pcorNetwork)

#' optimizedLambda returns the lambda value used in analysis.
#'
#' @param x A DNEAobject
#' @return The optimized lambda hyperparameter
#' @export
setGeneric("optimizedLambda",
           function(x) standardGeneric("optimizedLambda"))
setMethod("optimizedLambda", signature(x = "pcorNetwork"), function(x) x@hyperparameter$optimized_lambda)
#' #'
#' #'
#' #' @keywords internal
setGeneric("optimizedLambda<-",
           function(x,value) standardGeneric("optimizedLambda<-"))
setMethod("optimizedLambda<-",signature(x = "pcorNetwork"), function(x, value){

  x@hyperparameter$optimized_lambda <- value
  x
  validObject(x)
})

#' lambdas2Test returns the lambda value used in analysis.
#'
#' @param x A DNEAobject
#' @return The lambda values to evaluate in optimization
#' @keywords internal
lambdas2Test <- function(x){
  x@hyperparameter$tested_lambda_values
}
setGeneric("lambdas2Test",
           function(x) standardGeneric("lambdas2Test"))
setMethod("lambdas2Test", signature(x = "pcorNetwork"), function(x) x@hyperparameter$tested_lambda_values)
#' #'
#' #'
setGeneric("lambdas2Test<-",
           function(x,value) standardGeneric("lambdas2Test<-"))
setMethod("lambdas2Test<-",signature(x = "pcorNetwork"), function(x, value){

  x@hyperparameter$tested_lambda_values <- value
  x
  validObject(x)
})

#' BICscores returns the BIC scores for each lambda value evaluated
#'
#' @param x A DNEAobject
#' @return The optimized lambda hyperparameter
#' @keywords internal
BICscores.pcorNetwork <- function(x){
  x@hyperparameter$BIC_scores
}
setGeneric("BICscores",
           function(x) standardGeneric("BICscores"))
setMethod("BICscores", signature(x = "pcorNetwork"), function(x) BICscores.pcorNetwork)

#' includeMetadata adds info to metadata
#'
#' This function will take additional metadata and add it to the corresponding dataframe in the metadata
#' slot.
#'
#' @param object A DNEAobject
#' @param type sample or feature metadata
#' @param metadata a dataframe containing metadata to add
#'
#' @return A DNEAobject
#'
#' @export
includeMetadata <- function(object, type = c('sample', 'feature'), metadata){

  type = match.arg(type)
  if(type == 'sample'){
    if(all(sampleNames(object) == rownames(metadata))){
      for(i in 1:length(colnames(metadata))){


        object@metadata[["samples"]][[colnames(metadata)[i]]] <- metadata[, i]
      }
    } else{

      stop('new metadata order does not match sample order in DNEAobject')

    }
  } else{
    if(all(featureNames(object) == rownames(metadata)) |
       all(object@metadata[["features"]]$clean_feature_names == rownames(metadata))){
      for(i in 1:length(colnames(metadata))){

        new_metadata_colname <- colnames(metadata)[i]
        object@metadata[["features"]][[metadata[, i]]] <- metadata[, i]
      }
    } else{
      stop('new metadata order does not match feature order in DNEAobject')
    }
  }

  return(object)
}

#' getNetworkFiles will save the node and edge information
#'
#' This function will save the node and edge information as .csv files in the working directory.
#' The files are formatted for input into Cytoscape.
#'
#' @param object A DNEAobject
#'
#' @return The same DNEAobject as input and saves the node and edge information as .csv files in
#'          the working directory
#' @importFrom utils write.csv
#' @export
getNetworkFiles <- function(object, file_path){

  if(missing(file_path)){
    file_path <- paste0(getwd(), "/")
  }
  #save node list
  write.csv(object@node_list, paste0(file_path, object@project_name,'_nodelist_',Sys.Date(),'.csv'),
            row.names = FALSE)

  #save edge list
  write.csv(object@edge_list, paste0(file_path, object@project_name,'_edgelist_',Sys.Date(),'.csv'),
            row.names = FALSE)

  return(object)
}

