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
           netGSA_results = 'list'
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
    if(class(object@assays[[i]]) != 'matrix'){
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
      "  Project Name:  ", object@project_name, "\n",
      "  Un-scaled data: ", class(object@assays$Expression), "\n",
      "  Scaled data:  ", class(object@assays$NormalExpression), "\n",
      "  Samples:  ", paste0('There are ',length(object@metadata$Samples), ' samples.'), "\n",
      "  Features:  ", paste0('There are ',length(object@metadata$Features), 'Features.'), "\n",
      "  Conditions:  ", object@metadata$Conditions, "\n",
      "  Optimized Lambda: ", object@BIC[["optimizedLambda"]],
      "  Sub-clusters: ", unique(object@node_list$membership),
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
           function(x, mat) standardGeneric("expressionData<-"))
setMethod("expressionData<-", "DNEAobject", function(x, mat) {
  x@assays[['expresion_data']] <- mat
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
           function(x, mat) standardGeneric("scaledExpressionData<-"))
setMethod("scaledExpressionData<-", "DNEAobject", function(x, mat) {
  x@assays[['scaled_expression_data']] <- mat
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
           function(x,cond) standardGeneric("condition<-"))
setMethod("condition<-","DNEAobject", function(x,cond){
  x@metadata$condition_values <- cond
  validObject(x)
  return(x)
})

#'getter function for Samples
#'
#'@noRd
setGeneric("sampleNames",
           function(x) standardGeneric("sampleNames"))
setMethod("sampleNames", "DNEAobject", function(x) x@metadata$samples)

#'setter function for Samples
#'
#'@noRd
setGeneric("sampleNames<-",
           function(x,snames) standardGeneric("sampleNames<-"))
setMethod("sampleNames<-","DNEAobject", function(x,snames){
  x@metadata$samples <- snames
  validObject(x)
  return(x)
})

#'getter function for Features
#'
#'@noRd
setGeneric("featureNames",
           function(x) standardGeneric("featureNames"))
setMethod("featureNames", "DNEAobject", function(x) x@metadata$features)

#'getter function for clean Features
#'
#'@noRd
setGeneric("cleanFeatureNames",
           function(x) standardGeneric("cleanFeatureNames"))
setMethod("cleanFeatureNames", "DNEAobject", function(x) x@metadata$clean_feature_names)

#'setter function for condition
#'@import janitor
#'@noRd
setGeneric("featureNames<-",
           function(x,fnames) standardGeneric("featureNames<-"))
setMethod("featureNames<-","DNEAobject", function(x,fnames){
  x@metadata$features <- fnames
  x@metadata$clean_feature_names <- make_clean_names(fnames)
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




