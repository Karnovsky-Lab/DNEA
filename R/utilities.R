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
           Project.Name = 'character',
           Assays = 'list',
           Metadata = 'list',
           Dataset_summary = 'list',
           Nodes = 'data.frame',
           Edges = 'list',
           BIC = 'list',
           Adjacency.Matrices = 'list',
           Stable.Networks = 'list',
           Joint.Graph = 'igraph',
           NetGSA = 'list'
         )
)
#'Check Validity of "DNEAobject"
#'
#'@import methods
#'@noRd
setValidity("DNEAobject", function(object){
  if(length(object@Project.Name) > 1){
    "@Project.Name must be a character string"
  }
  for (i in length(object@Assays)){
    if(class(object@Assays[[i]]) != 'matrix'){
      "@Assays must be an expression matrix"
    }
    if(length(colnames(object@Assays[[i]])) != length(unique(colnames(object@Assays[[i]])))){
      "@Assays must be an expression matrix where each column is a unique feature."
    }
    if(!(is.numeric(object@Assays[[i]]))){
      "@Assays must be a matrix with numeric values."
    }
    if(all(rownames(object@Assays[[i]]) != object@Metadata$Samples)){
      "Samples are out of order"
    }
    if(all(colnames(object@Assays[[i]]) != object@Metadata$clean_Feature_Names)){
      "Features are out of order"
    }
  }
  if(class(object@Metadata$Samples) != 'character'){
    "@Metadata$Samples should be of class character"
  }
  if(class(object@Metadata$Features) != 'character'){
    " @Metadata$Features should be of class character"
  }
  if(class(object@Metadata$clean_Feature_Names) != 'character'){
    "@Metadata$clean_Feature_Names should be of class character"
  }
  if(class(object@Metadata$Condition) != 'factor' | length(levels(object@Metadata$Condition)) != 2){
    "@Metadata$Condition should be of class factor with 2 levels"
  }
})

#'modify show
#'
#'@noRd
setMethod("show", "DNEAobject", function(object) {
  cat(is(object)[[1]], "\n",
      "  Project Name:  ", object@Project.Name, "\n",
      "  Expression: ", class(object@Assays$Expression), "\n",
      "  NormalExpression:  ", class(object@Assays$NormalExpression), "\n",
      "  Sample:  ", paste0('There are ',length(object@Metadata$Samples), ' samples.'), "\n",
      "  Features:  ", paste0('There are ',length(object@Metadata$Features), 'Features.'), "\n",
      "  Conditions:  ", object@Metadata$Conditions, "\n",
      sep = ""
  )
})

#'getter function for Expression data
#'
#'@noRd
setGeneric("Expression",
           function(x) standardGeneric("Expression"))
setMethod("Expression", "DNEAobject", function(x) x@Assays$Expression)

#'setter function for Expression data
#'
#'@noRd
setGeneric("Expression<-",
           function(x, mat) standardGeneric("Expression<-"))
setMethod("Expression<-", "DNEAobject", function(x, mat) {
  x@Assays$Expression <- mat
  validObject(x)
  return(x)
})

#'getter function for NormalExpression data
#'
#'@noRd
setGeneric("NormalExpression",
           function(x) standardGeneric("NormalExpression"))
setMethod("NormalExpression", "DNEAobject", function(x) x@Assays$NormalExpression)

#'setter function for NormalExpression data
#'
#'@noRd
setGeneric("NormalExpression<-",
           function(x, mat) standardGeneric("NormalExpression<-"))
setMethod("NormalExpression<-", "DNEAobject", function(x, mat) {
  x@Assays$NormalExpression <- mat
  validObject(x)
  return(x)
})

#'getter function for Condition
#'
#'@noRd
setGeneric("condition",
           function(x) standardGeneric("condition"))
setMethod("condition", "DNEAobject", function(x) x@Metadata$Condition)

#'setter function for condition
#'
#'@noRd
setGeneric("condition<-",
           function(x,cond) standardGeneric("condition<-"))
setMethod("condition<-","DNEAobject", function(x,cond){
  x@Metadata$Condition <- cond
  validObject(x)
  return(x)
})

#'getter function for Samples
#'
#'@noRd
setGeneric("sampleNames",
           function(x) standardGeneric("sampleNames"))
setMethod("sampleNames", "DNEAobject", function(x) x@Metadata$Samples)

#'setter function for Samples
#'
#'@noRd
setGeneric("sampleNames<-",
           function(x,snames) standardGeneric("sampleNames<-"))
setMethod("sampleNames<-","DNEAobject", function(x,snames){
  x@Metadata$Samples <- snames
  validObject(x)
  return(x)
})

#'getter function for Features
#'
#'@noRd
setGeneric("featureNames",
           function(x) standardGeneric("featureNames"))
setMethod("featureNames", "DNEAobject", function(x) x@Metadata$Features)

#'getter function for clean Features
#'
#'@noRd
setGeneric("cleanFeatureNames",
           function(x) standardGeneric("cleanFeatureNames"))
setMethod("cleanFeatureNames", "DNEAobject", function(x) x@Metadata$clean_Feature_Names)

#'setter function for condition
#'@import janitor
#'@noRd
setGeneric("featureNames<-",
           function(x,fnames) standardGeneric("featureNames<-"))
setMethod("featureNames<-","DNEAobject", function(x,fnames){
  x@Metadata$Features <- fnames
  x@Metadata$clean_Feature_Names <- make_clean_names(fnames)
  validObject(x)
  return(x)
})

#'getter function for number of features
#'
#'@noRd
setGeneric("numFeatures",
           function(x) standardGeneric("numFeatures"))
setMethod("numFeatures", "DNEAobject", function(x) x@Dataset_summary$num_features)

#'getter function for number of samples
#'
#'@noRd
setGeneric("numSamples",
           function(x) standardGeneric("numSamples"))
setMethod("numSamples", "DNEAobject", function(x) x@Dataset_summary$num_samples)


