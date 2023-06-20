

#' @include all-classes.R
#' @rdname expressionData
#' @export
setGeneric("expressionData", function(x, type) standardGeneric("expressionData"))


#' @rdname conditionLevels
#' @export
setGeneric("conditionLevels", function(x) standardGeneric("conditionLevels"))


#' @rdname conditions
#' @export
setGeneric("conditions", function(x) standardGeneric("conditions"))


#' @rdname conditions
#' @export
setGeneric("conditions<-", function(x, value) standardGeneric("conditions<-"))


#' @rdname sampleNames
#' @export
setGeneric("sampleNames", function(x) standardGeneric("sampleNames"))


#' @rdname featureNames
#' @export
setGeneric("featureNames", function(x, original = TRUE) standardGeneric("featureNames"))


#' @rdname numFeatures
#' @export
setGeneric("numFeatures", function(x) standardGeneric("numFeatures"))


#' @rdname numSamples
#' @export
setGeneric("numSamples", function(x) standardGeneric("numSamples"))


#' @rdname optimizedLambda
#' @export
setGeneric("optimizedLambda", function(x) standardGeneric("optimizedLambda"))


#' @rdname optimizedLambda
#' @export
setGeneric("optimizedLambda<-", function(x, value) standardGeneric("optimizedLambda<-"))


#' @rdname lambdas2Test
#' @export
setGeneric("lambdas2Test", function(x) standardGeneric("lambdas2Test"))


#' @rdname lambdas2Test
#' @export
setGeneric("lambdas2Test<-", function(x, value) standardGeneric("lambdas2Test<-"))


#' @rdname BICscores
#' @export
setGeneric("BICscores", function(x) standardGeneric("BICscores"))










