
#' set "igraph" class for consensus cluster
#'
#' @import igraph
#' @keywords internal
#' @noRd
setOldClass('igraph')

#'Set "DNEAobject" class
#'
#' @import methods
#' @noRd
setClass(Class = "DNEAobject",
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
           joint_graph = 'igraph',
           netGSA_results = 'list',
           feature_membership = 'list'
         )
)
#'Check Validity of "DNEAobject"
#'
#' @import methods
#' @noRd
setValidity("DNEAobject", function(object){
  if(length(object@project_name) > 1){
    "@project_name must be a character string"
  }
  for (i in length(object@assays)){
    if(class(object@assays[[i]])[[1]] != 'matrix'){
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
  if(class(object@metadata$samples) != 'character'){
    "@metadata$Samples should be of class character"
  }
  if(class(object@metadata$features) != 'character'){
    " @metadata$Features should be of class character"
  }
  if(class(object@metadata$clean_feature_names) != 'character'){
    "@metadata$clean_Feature_Names should be of class character"
  }
  if(class(object@metadata$condition_values) != 'factor' | length(levels(object@metadata$condition_values)) != 2){
    "@metadata$Condition should be of class factor with 2 levels"
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
#' @return the expression matrix
#' @export
setGeneric("expressionData",
           function(x) standardGeneric("expressionData"))
setMethod("expressionData", "DNEAobject", function(x) x@assays[['expression_data']])


# setGeneric("expressionData<-",
#            function(x, value) standardGeneric("expressionData<-"))
# setMethod("expressionData<-", "DNEAobject", function(x, value) {
#   x@assays[['expresion_data']] <- value
#   validObject(x)
#   return(x)
# })

#' scaledExpressionData retrieves the scaled expression data from the assays slot.
#'
#' This function takes in a DNEAobject and return the scaled expression matrix located in the
#' assays slot of the object. You can also use it to input scaled expression data into the DNEAobject
#'
#' @param x A DNEAobject
#' @return the scaled expression matrix
#' @export
setGeneric("scaledExpressionData",
           function(x) standardGeneric("scaledExpressionData"))
setMethod("scaledExpressionData", "DNEAobject", function(x) x@assays[['scaled_expression_data']])

##'setter function for NormalExpression data
#
#setGeneric("scaledExpressionData<-",
#           function(x, value) standardGeneric("scaledExpressionData<-"))
#setMethod("scaledExpressionData<-", "DNEAobject", function(x, value) {
#  x@assays[['scaled_expression_data']] <- value
#  validObject(x)
#  return(x)
#})

#' condition retrieves the condition values for each sample from the metadata slot.
#'
#' This function takes in a DNEAobject and return the condition values located in the
#' metatdata slot of the object. You can also use it to input conditions data into the DNEAobject
#'
#' @param x A DNEAobject
#' @return A vector of the condition values
#' @export
setGeneric("condition",
           function(x) standardGeneric("condition"))
setMethod("condition", "DNEAobject", function(x) x@metadata$condition_values)

##'setter function for condition

#setGeneric("condition<-",
#           function(x,value) standardGeneric("condition<-"))
#etMethod("condition<-","DNEAobject", function(x,value){
#  x@metadata$condition_values <- value
#  validObject(x)
#  return(x)
#})

#' sampleNames retrieves the sample names from the metadata slot.
#'
#' This function takes in a DNEAobject and return the sample names located in the
#' metadata slot of the object.
#'
#' @param x A DNEAobject
#' @return A vector of sample names
#' @export
setGeneric("sampleNames",
           function(x) standardGeneric("sampleNames"))
setMethod("sampleNames", "DNEAobject", function(x) x@metadata[["samples"]]$samples)
#
##'setter function for Samples
##'
#setGeneric("sampleNames<-",
#           function(x,value) standardGeneric("sampleNames<-"))
#setMethod("sampleNames<-","DNEAobject", function(x,value){
#  x@metadata[["samples"]]$samples <- value
#  x@assays[[""]]
#  validObject(x)
#  return(x)
#})

#' featureNames retrieves the feature names from the metadata slot.
#'
#' This function takes in a DNEAobject and return the feature names located in the
#' metadata slot of the object.
#'
#' @param x A DNEAobject
#' @return A vector of feature names
#' @export
setGeneric("featureNames",
           function(x) standardGeneric("featureNames"))
setMethod("featureNames", "DNEAobject", function(x) x@metadata[["features"]]$features)

#' cleanFeatureNames retrieves the feature names from clean_feature names in the metadata slot.
#'
#' This function takes in a DNEAobject and return the feature names modified to avoid errors with
#' native R clashes. The names are located in the metadata slot of the object.
#'
#' @param x A DNEAobject
#' @return A vector of modified feature names
#' @export
setGeneric("cleanFeatureNames",
           function(x) standardGeneric("cleanFeatureNames"))
setMethod("cleanFeatureNames", "DNEAobject", function(x) x@metadata[["features"]]$clean_feature_names)


# setGeneric("featureNames<-",
#           function(x,value) standardGeneric("featureNames<-"))
# setMethod("featureNames<-","DNEAobject", function(x,value){
#  x@metadata[["features"]]$features <- value
#  x@metadata[["features"]]$clean_feature_names <- make_clean_names(value)
#  validObject(x)
#  return(x)
# })

#' numFeatures finds the total number of features in the dataset
#'
#' @param x A DNEAobject
#' @return The number of features in the dataset
#' @export
setGeneric("numFeatures",
           function(x) standardGeneric("numFeatures"))
setMethod("numFeatures", "DNEAobject", function(x) x@dataset_summary$num_features)

#' numSamples finds the total number of samples in the dataset
#'
#' @param x A DNEAobject
#' @return The number of samples in the dataset
#' @export
setGeneric("numSamples",
           function(x) standardGeneric("numSamples"))
setMethod("numSamples", "DNEAobject", function(x) x@dataset_summary$num_samples)

#' optimizedLambda returns the lambda value used in analysis.
#'
#' @param x A DNEAobject
#' @return The optimized lambda hyperparameter
#' @export
setGeneric("optimizedLambda",
           function(x) standardGeneric("optimizedLambda"))
setMethod("optimizedLambda", "DNEAobject", function(x) x@hyperparameter$optimized_lambda)


setGeneric("optimizedLambda<-",
           function(x,value) standardGeneric("optimizedLambda<-"))
setMethod("optimizedLambda<-","DNEAobject", function(x, value){

  x@hyperparameter$optimized_lambda <- value
  x
  validObject(x)
})

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

        new_metadata_colname <- colnames(metadata)[i]
        object@metadata[["samples"]] <- data.frame(object@metadata[["samples"]],
                                                   new_metadata_colname = metadata[,i])
      }
    } else{

      stop('new metadata order does not match sample order in DNEAobject')

    }
  } else{
    if(all(featureNames(object) == rownames(metadata)) |
       all(object@metadata[["features"]]$clean_feature_names == rownames(metadata))){
      for(i in 1:length(colnames(metadata))){

        new_metadata_colname <- colnames(metadata)[i]
        object@metadata[["features"]] <- data.frame(object@metadata[["features"]],
                                                    new_metadata_colname = metadata[,i])
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
#'  @param object A DNEAobject
#'
#'  @return The same DNEAobject as input and saves the node and edge information as .csv files in
#'          the working directory
#' @importFrom utils write.csv
#' @export
getNetworkFiles <- function(object){

  #save node list
  write.csv(object@node_list, paste0(object@project_name,'_nodelist_',Sys.Date(),'.csv'), row.names = FALSE)

  #save edge list
  write.csv(object@edge_list, paste0(object@project_name,'_edgelist_',Sys.Date(),'.csv'), row.names = FALSE)

  return(object)
}

