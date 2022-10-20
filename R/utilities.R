#' set "igraph" class for consensus cluster
#'
#' @import igraph
#'@noRd
setOldClass('igraph')

#'Set "DNEAobject" class
#'
#'@import methods
#'@noRd
#Internal Function to create DNEAobject
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
#'@import methods
#'@noRd
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

#'modify show
#'
#'@noRd
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

#'split data by condition
#'
#'@noRd
split_by_condition <- function(dat, condition_levels, condition_by_sample){

  #create key for separating the data by key and running diagnostic tests, feature DE calculations
  separated_conditions_data <- vector(mode = 'list', length = length(condition_levels))
  names(separated_conditions_data) <- condition_levels

    for(cond in condition_levels){
      separated_conditions_data[[cond]] <- t(dat)[,(condition_by_sample == cond)]

    }
  return(separated_conditions_data)

}
#'getter function for Expression data
#'
#'@noRd
setGeneric("expressionData",
           function(x) standardGeneric("expressionData"))
setMethod("expressionData", "DNEAobject", function(x) x@assays[['expression_data']])

#'setter function for Expression data
#'
#'@noRd
setGeneric("expressionData<-",
           function(x, value) standardGeneric("expressionData<-"))
setMethod("expressionData<-", "DNEAobject", function(x, value) {
  x@assays[['expresion_data']] <- value
  validObject(x)
  return(x)
})

#'getter function for NormalExpression data
#'
#'@noRd
setGeneric("scaledExpressionData",
           function(x) standardGeneric("scaledExpressionData"))
setMethod("scaledExpressionData", "DNEAobject", function(x) x@assays[['scaled_expression_data']])

#'setter function for NormalExpression data
#'
#'@noRd
setGeneric("scaledExpressionData<-",
           function(x, value) standardGeneric("scaledExpressionData<-"))
setMethod("scaledExpressionData<-", "DNEAobject", function(x, value) {
  x@assays[['scaled_expression_data']] <- value
  validObject(x)
  return(x)
})

#'getter function for Condition
#'
#'@noRd
setGeneric("condition",
           function(x) standardGeneric("condition"))
setMethod("condition", "DNEAobject", function(x) x@metadata$condition_values)

#'setter function for condition
#'
#'@noRd
setGeneric("condition<-",
           function(x,value) standardGeneric("condition<-"))
setMethod("condition<-","DNEAobject", function(x,value){
  x@metadata$condition_values <- value
  validObject(x)
  return(x)
})

#'getter function for Samples
#'
#'@noRd
setGeneric("sampleNames",
           function(x) standardGeneric("sampleNames"))
setMethod("sampleNames", "DNEAobject", function(x) x@metadata[["samples"]]$samples)

#'setter function for Samples
#'
#'@noRd
setGeneric("sampleNames<-",
           function(x,value) standardGeneric("sampleNames<-"))
setMethod("sampleNames<-","DNEAobject", function(x,value){
  x@metadata[["samples"]]$samples <- value
  validObject(x)
  return(x)
})

#'getter function for Features
#'
#'@noRd
setGeneric("featureNames",
           function(x) standardGeneric("featureNames"))
setMethod("featureNames", "DNEAobject", function(x) x@metadata[["features"]]$features)

#'getter function for clean Features
#'
#'@noRd
setGeneric("cleanFeatureNames",
           function(x) standardGeneric("cleanFeatureNames"))
setMethod("cleanFeatureNames", "DNEAobject", function(x) x@metadata[["features"]]$clean_feature_names)

#'setter function for condition
#'#'@importFrom janitor make_clean_names
#'@noRd
setGeneric("featureNames<-",
           function(x,value) standardGeneric("featureNames<-"))
setMethod("featureNames<-","DNEAobject", function(x,value){
  x@metadata[["features"]]$features <- value
  x@metadata[["features"]]$clean_feature_names <- make_clean_names(value)
  validObject(x)
  return(x)
})

#'getter function for number of features
#'
#'@noRd
setGeneric("numFeatures",
           function(x) standardGeneric("numFeatures"))
setMethod("numFeatures", "DNEAobject", function(x) x@dataset_summary$num_features)

#'getter function for number of samples
#'
#'@noRd
setGeneric("numSamples",
           function(x) standardGeneric("numSamples"))
setMethod("numSamples", "DNEAobject", function(x) x@dataset_summary$num_samples)

#'getter function for optimized lambda value
#'
#'@noRd
setGeneric("optimizedLambda",
           function(x) standardGeneric("optimizedLambda"))
setMethod("optimizedLambda", "DNEAobject", function(x) x@hyperparameter$optimized_lambda)

#'setter function for optimized lambda value
#'
#'@noRd
setGeneric("optimizedLambda<-",
           function(x,value) standardGeneric("optimizedLambda<-"))
setMethod("optimizedLambda<-","DNEAobject", function(x, value){

  x@hyperparameter$optimized_lambda <- value
  x
  validObject(x)
})

#add info to metadata
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


