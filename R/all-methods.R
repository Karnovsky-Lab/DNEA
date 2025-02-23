#' @include all-classes.R
#' @include all-generics.R
NULL

#' Retrieve information about a DNEA object
#'
#' @describeIn DNEA-class
#' This function will display a summary of the information
#' stored within a \code{\link{DNEA}} object.
#'
#'
#' @param object A \code{\link{DNEA}} object.
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}},
#' \code{\link{aggregateFeatures}}
#'
#' @returns A summary of the information stored in a
#' \code{\link{DNEA}} object.
#'
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEA object
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' dnw
#' @export
setMethod("show", "DNEA", function(object) {
  cat(is(object)[[1]],
      "\n  Project Name -  ", object@project_name,
      "\n  Samples -  ", 'There are ',length(sampleNames(object)),
      ' samples.',
      "\n  Features -  ", 'There are ',length(featureNames(object)),
      ' Features.',
      "\n  Conditions -  ", "control: ", networkGroups(object)[[1]],
                            " case: ", networkGroups(object)[[2]],
      "\n  Optimized Lambda - ", ifelse(is.null(optimizedLambda(object)),
                                        "Parameter tuning not performed",
                                        paste0("lambda: ",
                                               optimizedLambda(object),
                                               " will be used in analysis.")),
      "\n  Metabolic modules: ", ifelse(is.null(nodeList(object)[["membership"]]),
                                        "Consensus Clustering not performed",
                                        paste0("there are ",
                                               length(unique(nodeList(object)[["membership"]])),
                                               " subnetworks.")),
      "\n  netGSA results: ", ifelse(nrow(netGSAresults(object)) == 0,
                                     "Enrichment analysis not performed",
                                     paste0(sum(netGSAresults(object)$NetGSA_pFDR <= 0.05),
                                            " enriched modules\n")),
      sep=""
  )
})

#' Retrieve information about a DNEAinputSummary object
#'
#' @describeIn DNEAinputSummary-class
#' This function will display the number of samples, number of features,
#' and diagnostics values of the input data set to a
#' \code{\link{DNEA}} object.
#'
#' @param object A DNEAinputSummary object
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}},
#' \code{\link{aggregateFeatures}}
#'
#' @returns A summary of the input data to \code{\link{createDNEAobject}}.
#'
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEA object
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' datasetSummary(dnw)
#' @export
setMethod("show", "DNEAinputSummary", function(object){
  cat(is(object)[[1]],
             "\n  Number of Samples  -  ", numSamples(object),
             "\n  Number of Features  -  ", numFeatures(object),"\n",
      sep="")
  print(as.matrix(diagnostics(object)))
})

projectName.DNEA <- function(x){

  x@project_name
}
#' Return the name of the current experiment
#'
#' This function returns the name of the DNEA experiment.
#'
#' @param x A \code{\link{DNEA}} object.
#'
#' @returns The name of the DNEA experiment.
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEA object
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' projectName(dnw)
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}}
#' @rdname projectName-methods
#' @aliases projectName
#' @export
setMethod("projectName", signature(x= "DNEA"),
          projectName.DNEA)

expressionData.DNEA <- function(x, assay=names(assays(x))){
  assay <- match.arg(assay)
  output <- assays(x)[[assay]]
  return(output)
}
#' Access expression data within a DNEA object,
#'
#' This function accesses the expression data stored in the
#' assays slot of the \code{\link{DNEA}} object. The
#' output is an \emph{n x m} matrix with one row for each
#' sample and one column for each feature in the data.
#'
#'
#' @param x A \code{\link{DNEA}} object.
#'
#' @param assay A character string corresponding to the
#' data to retrieve: "input_data" retrieves the data as
#' it was input, "log_input_data" retrieves the input data
#' after log transforming, and "log-scaled_data"
#' retrieves a list of matrices corresponding to the
#' log-scaled data for each experimental condition,
#' respectively. Any other externally transformed
#' data that is stored in the DNEA object can be accessed
#' by providing its name to the assay parameter.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}},\code{\link{aggregateFeatures}}

#' @returns The expression matrix specified by the user.
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEA object
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' expressionData(x=dnw, assay="input_data")
#' expressionData(x=dnw, assay="log_input_data")
#' expressionData(x=dnw, assay="log-scaled_data")
#' @rdname expressionData-methods
#' @aliases expressionData
#' @export
setMethod("expressionData",signature(x="DNEA"),
          expressionData.DNEA)

assays.DNEA <- function(x){

  x@assays
}
#' @rdname assays-methods
#' @aliases assays
#' @keywords internal
#' @noRd
setMethod("assays", signature(x="DNEA"),
          assays.DNEA)

assaysReplace.DNEA <- function(x, value){
  x@assays <- value
  validObject(x)
  x
}
#' @rdname assays-methods
#' @aliases assays
#' @keywords internal
#' @noRd
setReplaceMethod("assays", signature(x="DNEA"),
                 assaysReplace.DNEA)
networkGroupIDs.DNEA <- function(x){

  x@metadata$network_group_IDs
}

#' Access and set the experimental group labels
#'
#' This function accesses the experimental group labels for
#' each sample stored in the metadata slot of a
#' \code{\link{DNEA}} object.
#'
#' @param x A \code{\link{DNEA}} object.
#'
#' @param value a character string name corresponding to a column
#' name of the sample metadata data frame.
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{includeMetadata}}
#' @returns A vector of the unique condition labels.
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEA object
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' networkGroupIDs(dnw)
#' @rdname networkGroupIDs-methods
#' @aliases networkGroupIDs
#' @export
setMethod("networkGroupIDs", signature(x="DNEA"),
          networkGroupIDs.DNEA)

networkGroups.DNEA <- function(x){

  x@metadata$network_groups
}
#' Retrieve the unique group values of the experimental condition
#'
#' This function takes in a \code{\link{DNEA}} object and
#' returns the unique group labels of the experimental
#' condition in the data set.
#'
#' @param x A \code{\link{DNEA}} or
#' \code{\link{collapsed_DNEA}}
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{networkGroupIDs}},
#' \code{\link{createDNEAobject}}
#' @returns A vector of the condition values.
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEA object
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' networkGroups(dnw)
#' @rdname networkGroups-methods
#' @aliases networkGroups
#' @export
setMethod("networkGroups", signature(x="DNEA"),
          networkGroups.DNEA)

metaData.DNEA <- function(x, type=c("samples", "features")){

  type <- match.arg(type)
  if(type == "samples"){
    x@metadata[["samples"]]
  }else if(type == "features"){
    x@metadata[["features"]]
  }
}
#' Retrieve metadata stored in a DNEA
#'
#' This function retrieves the specified metadata stored
#' in the metadata slot of the
#' \code{\link{DNEA}} object.
#'
#' @param x A \code{\link{DNEA}} object.
#' @param type A character string indicating the type
#' of metadata to access. Can be "sample" or "feature".
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}}, \code{\link{includeMetadata}}
#' @returns A data frame of the indicated metadata
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEA
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' metaData(dnw, type = "sample")
#' @seealso \code{\link{includeMetadata}}
#' @rdname metaData-methods
#' @aliases metaData
#' @export
setMethod("metaData", signature(x="DNEA"),
          metaData.DNEA)

metaDataReplace.DNEA <- function(x,
                                    type=c("samples", "features"),
                                    value){

  type <- match.arg(type)
  x@metadata[[type]] <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("metaData", signature(x="DNEA"),
                 metaDataReplace.DNEA)

sampleNames.DNEA <- function(x){

  x@metadata[["samples"]]$samples
}
#' Retrieve the sample names from the metadata slot.
#'
#' This function accesses the sample names stored in the
#' metadata slot of the \code{\link{DNEA}} object.
#'
#' @param x A \code{\link{DNEA}} object.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}}
#' @returns A character vector of sample names.
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEA object
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' sampleNames(dnw)
#' @rdname sampleNames-methods
#' @aliases sampleNames
#' @export
setMethod("sampleNames", signature(x="DNEA"),
          sampleNames.DNEA)

featureNames.DNEA <- function(x, original=FALSE){

  if(isTRUE(original)){
    x@metadata[["features"]]$feature_names
  }else if(isFALSE(original)){
    x@metadata[["features"]]$clean_feature_names
  }
}
#' Retrieve the feature names from the metadata slot.
#'
#' This function accesses the feature names stored in the
#' metadata slot of the \code{\link{DNEA}} object.
#'
#' @param x A \code{\link{DNEA}} object.
#' @param original "TRUE" returns the original feature names
#' and "FALSE" returns the feature names that have been
#' modified to avoid errors as a result of special characters
#' using \code{\link[janitor:make_clean_names]{make_clean_names}}.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}}
#' @returns A character vector of feature names.
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEA object
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' featureNames(dnw, original=TRUE)
#' @rdname featureNames-methods
#' @aliases featureNames
#' @export
setMethod("featureNames", signature(x="DNEA"),
          featureNames.DNEA)

numFeatures.DNEA <- function(x){
  x@dataset_summary@num_features
}
#' Retrieve the total number of features in the dataset
#'
#' This function prints to console the total number of
#' features in the data set
#'
#' @param x A \code{\link{DNEA}} object.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}}
#' @returns The number of features in the data set.
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEA object
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' numFeatures(dnw)
#' @rdname numFeatures-methods
#' @aliases numFeatures
#' @export
setMethod("numFeatures", signature(x="DNEA"),
          numFeatures.DNEA)

numFeatures.DNEAinputSummary <- function(x){

  x@num_features
}
#' @rdname numFeatures-methods
#' @export
setMethod("numFeatures", signature(x="DNEAinputSummary"),
          numFeatures.DNEAinputSummary)

numSamples.DNEA <- function(x){

  x@dataset_summary@num_samples
}
#' Retrieves the total number of samples in the dataset
#'
#' This function prints to console the total number of
#' samples in the data set.
#'
#' @param x A \code{\link{DNEA}} object.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}}
#' @returns The number of samples in the data set.
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEA object
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' numSamples(dnw)
#' @rdname numSamples-methods
#' @aliases numSamples
#' @export
setMethod("numSamples", signature(x="DNEA"),
          numSamples.DNEA)

numSamples.DNEAinputSummary <- function(x){

  x@num_samples
}
#' @rdname numSamples-methods
#' @aliases numSamples
#' @export
setMethod("numSamples", signature(x="DNEAinputSummary"),
          numSamples.DNEAinputSummary)

optimizedLambda.DNEA <- function(x){

  x@hyperparameter$optimized_lambda
}
#' Access the lambda value used in analysis
#'
#' The function takes as input a \code{\link{DNEA}} object and returns
#' the hyper parameter (lambda) that is currently being used for the analysis.
#' The user may also provide a single-value numeric vector to change the
#' lambda value for analysis.
#'
#' @param x A \code{\link{DNEA}} object.
#'
#' @param value a single-value numeric vector corresponding to the lambda
#' value to use in analysis.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{BICtune}}
#' @returns The optimized lambda hyperparameter.
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEA object
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' optimizedLambda(dnw) <- 0.15
#' optimizedLambda(dnw)
#' @rdname optimizedLambda-methods
#' @aliases optimizedLambda
#' @export
setMethod("optimizedLambda", signature(x="DNEA"),
          optimizedLambda.DNEA)

optimizedLambdaReplace.DNEA <- function(x, value){

  x@hyperparameter$optimized_lambda <- value
  validObject(x)
  x
}
#' @rdname optimizedLambda-methods
#' @aliases optimizedLambda
#' @export
setReplaceMethod("optimizedLambda", signature(x="DNEA"),
                 optimizedLambdaReplace.DNEA)

lambdas2Test.DNEA <- function(x){

  x@hyperparameter$tested_lambda_values
}
#' Access the lambda values tested during
#' hyper parameter optimization
#'
#' The function takes as input a \code{\link{DNEA}}
#' object and returns the lambda values that were testing
#' during hyper parameter optimization performed via
#' \code{\link{BICtune}}.
#'
#'
#' @param x A \code{\link{DNEA}} or
#' \code{\link{collapsed_DNEA}} object.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{BICtune}}
#' @returns The lambda values to evaluate in optimization.
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEA object
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' lambdas2Test(dnw)
#' @rdname lambdas2Test-methods
#' @aliases lambdas2Test
#' @export
setMethod("lambdas2Test", signature(x="DNEA"),
          lambdas2Test.DNEA)

lambdas2TestReplace.DNEA <- function(x, value){

  x@hyperparameter$tested_lambda_values <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("lambdas2Test", signature(x="DNEA"),
                 lambdas2TestReplace.DNEA)

BICscores.DNEA <- function(x){

  x@hyperparameter$BIC_scores
}
#' Access the BIC scores for each lambda value evaluated
#'
#' The function takes as input a \code{\link{DNEA}}
#' object and returns the BIC values for each lambda tested
#' during hyper parameter optimization performed via
#' \code{\link{BICtune}}.
#'
#'
#' @param x A \code{\link{DNEA}} object.
#' @param value a list of two lists that consist of the
#' likelihood and BIC scores for each tested lambda value.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{BICtune}}
#' @returns The optimized lambda hyperparameter.
#' @examples
#' #dnw is a DNEA with the results generated for the example data
#' #accessed by running data(TEDDY) in the console. The workflow
#' #for this data can be found in the vignette accessed by
#' #running browseVignettes("DNEA") in the console.
#' data(dnw)
#'
#' BICscores(dnw)
#' @rdname BICscores-methods
#' @aliases BICscores
#' @export
setMethod("BICscores", signature(x="DNEA"),
          BICscores.DNEA)

BICscoresReplace.DNEA <- function(x, value){

  x@hyperparameter$BIC_scores <- value
  validObject(x)
  x
}
#' @rdname BICscores-methods
#' @aliases BICscores
setReplaceMethod("BICscores", signature(x="DNEA"),
                 BICscoresReplace.DNEA)

selectionResults.DNEA <- function(x){

  x@stable_networks$selection_results
}

#' Access and set the edge selection results from stabilitySelection()
#'
#' The function takes as input a \code{\link{DNEA}} object and returns
#' an \emph{m x n} matrix of selection results for every possible network
#' edge calculated during \code{\link{stabilitySelection}}.
#'
#'
#' @param x A \code{\link{DNEA}} object.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{stabilitySelection}},\code{\link{selectionProbabilities}}
#'
#' @returns A \code{\link{DNEA}} object after filling the
#' selection_results section of the stable_networks slot.
#'
#' @examples
#' #dnw is a DNEA with the results generated for the example data
#' #accessed by running data(TEDDY) in the console. The workflow
#' #for this data can be found in the vignette accessed by
#' #running browseVignettes("DNEA") in the console.
#' data(dnw)
#'
#' selectionResults(dnw)
#' @rdname selectionResults-methods
#' @aliases selectionResults
#' @export
setMethod("selectionResults", signature(x="DNEA"),
          selectionResults.DNEA)

selectionResultsReplace.DNEA <- function(x, value){

  x@stable_networks$selection_results <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("selectionResults", signature(x="DNEA"),
                 selectionResultsReplace.DNEA)

selectionProbabilities.DNEA <- function(x){

  x@stable_networks$selection_probabilities
}

#' Access and set the edge selection probabilities from stabilitySelection()
#'
#' The function takes as input a \code{\link{DNEA}} object and returns an
#' \emph{m x n} matrix of selection probabilities for every possible network
#' edge calculated during \code{\link{stabilitySelection}}.
#'
#'
#' @param x A \code{\link{DNEA}} object.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{stabilitySelection}},\code{\link{selectionResults}}
#'
#' @returns A \code{\link{DNEA}} object after filling the
#' selection_probabilities section of the stable_networks slot.
#' @examples
#' #dnw is a \code{\link{DNEA}} object with the results generated
#' #for the example data accessed by running data(TEDDY) in the
#' #console. The workflow for this data can be found in the
#' #vignette accessed by running browseVignettes("DNEA")
#' #in the console.
#' data(dnw)
#'
#' selectionProbabilities(dnw)
#' @rdname selectionProbabilities-methods
#' @aliases selectionProbabilities
#' @export
setMethod("selectionProbabilities", signature(x="DNEA"),
          selectionProbabilities.DNEA)

selectionProbabilitiesReplace.DNEA <- function(x, value){

  x@stable_networks$selection_probabilities <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("selectionProbabilities", signature(x="DNEA"),
                 selectionProbabilitiesReplace.DNEA)

edgeList.DNEA <- function(x){

  x@edge_list
}
#' Access the edge list
#'
#' The function takes as input a \code{\link{DNEA}} object and
#' returns the edge list created by the \code{\link{getNetworks}}
#' function.
#'
#' @param x a \code{\link{DNEA}} object.
#'
#' @param value a data frame of edges in the network.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{getNetworks}},\code{\link{filterNetworks}},
#' \code{\link{getNetworkFiles}}
#'
#' @returns A data frame corresponding to the edge list
#' determined by DNEA.
#' @examples
#' #dnw is a \code{\link{DNEA}} object with the results
#' #generated for the example data accessed by running
#' #data(TEDDY) in the console. The workflow for this data
#' #can be found in the vignette accessed by running
#' #browseVignettes("DNEA") in the console.
#' data(dnw)
#'
#' edgeList(dnw)
#' @rdname edgeList-methods
#' @aliases edgeList
#' @export
setMethod("edgeList", signature(x="DNEA"),
          edgeList.DNEA)

edgeListReplace.DNEA <- function(x, value){

  x@edge_list <- value
  validObject(x)
  x
}
#' @rdname edgeList-methods
#' @aliases edgeList
setReplaceMethod("edgeList", signature(x="DNEA"),
                 edgeListReplace.DNEA)

nodeList.DNEA <- function(x){

  x@node_list
}
#' Access the node list
#'
#' The function takes as input a \code{\link{DNEA}} object
#' and returns the node list created from the
#' \code{\link{createDNEAobject}} function.
#'
#' @param x a \code{\link{DNEA}} object.
#'
#' @param value a data frame of nodes in the network.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}},\code{\link{clusterNet}},
#' \code{\link{getNetworkFiles}}
#'
#' @returns A data frame corresponding to the node
#' list determined by DNEA.
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEA object
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' nodeList(dnw)
#' @rdname nodeList-methods
#' @aliases nodeList
#' @export
setMethod("nodeList", signature(x="DNEA"),
          nodeList.DNEA)

nodeListReplace.DNEA <- function(x, value){

  x@node_list <- value
  validObject(x)
  x
}
#' @rdname nodeList-methods
#' @aliases nodeList
setReplaceMethod("nodeList", signature(x="DNEA"),
                 nodeListReplace.DNEA)

diagnostics.DNEA <- function(x){

  x@dataset_summary@diagnostic_values
}
#' Retrieve the diagnostic values for the input expression data
#'
#' This function retrieves the diagnostic values calculated
#' for the input expression data by the
#' \code{\link{createDNEAobject}} function.
#'
#' @param x \code{\link{DNEA}} or
#' \code{DNEAinputSummary} object.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}}, \code{\link{aggregateFeatures}}
#'
#' @returns Returns the diagnostic values for
#' the input expression data.
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEA object
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' diagnostics(dnw)
#' @rdname diagnostics-methods
#' @aliases diagnostics
#' @export
setMethod("diagnostics", signature(x="DNEA"),
          diagnostics.DNEA)

diagnosticsReplace.DNEA <- function(x, value){

  x@dataset_summary$diagnostic_values <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("diagnostics", signature(x="DNEA"),
                 diagnosticsReplace.DNEA)

diagnostics.DNEAinputSummary <- function(x){

  x@diagnostic_values
}
#' @rdname diagnostics-methods
#' @aliases diagnostics
setMethod("diagnostics", signature(x="DNEAinputSummary"),
          diagnostics.DNEAinputSummary)

diagnosticsReplace.DNEAinputSummary <- function(x, value){

  x@diagnostic_values <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("diagnostics", signature(x="DNEAinputSummary"),
                 diagnosticsReplace.DNEAinputSummary)

datasetSummary.DNEA <- function(x){

  x@dataset_summary
}
#' Access the dataset_summary slot of a DNEA object
#'
#' This function prints to console the number of samples, number of
#' features, and diagnostic values of the input data stored in the
#' dataset_summary slot of the \code{\link{DNEA}}.
#'
#'
#' @param x A \code{\link{DNEA}} object.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{createDNEAobject}},\code{\link{aggregateFeatures}}
#'
#' @returns The numbers of samples/features and diagnostic values
#' of the input data stored in the dataset_summary slot of the
#' \code{\link{DNEA}}.
#' @examples
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEA object
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' datasetSummary(dnw)
#' @rdname datasetSummary-methods
#' @aliases datasetSummary
#' @export
setMethod("datasetSummary", signature(x="DNEA"),
          datasetSummary.DNEA)

datasetSummaryReplace.DNEA <- function(x, value){

  x@dataset_summary <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("datasetSummary", signature(x="DNEA"),
                 datasetSummaryReplace.DNEA)

adjacencyMatrix.DNEA <- function(x, weighted=FALSE){

  if(weighted){
    x@adjacency_matrix$weighted_adjacency
  } else if(!weighted){
    x@adjacency_matrix$unweighted_adjacency
  }
}
#' Retrieve the weighted or unweighted adjacency matrix
#'
#' The function takes as input a \code{\link{DNEA}} object
#' and returns the weighted or un-weighted adjacency matrix for
#' each group network constructed via the
#' \code{\link{getNetworks}} function.
#'
#'
#' @param x A \code{\link{DNEA}} object.
#' @param weighted TRUE/FALSE indicating whether the
#' weighted unweighted adjacency matrix should be
#' returned.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{getNetworks}}
#' @returns A matrix corresponding to the adjacency
#' matrix specified.
#' @examples
#' #dnw is a \code{\link{DNEA}} object with the results
#' #generated for the example data accessed by running
#' #data(TEDDY) in the console. The workflow for this data
#' #can be found in the vignette accessed by running
#' #browseVignettes("DNEA") in the console.
#'
#' adjacencyMatrix(dnw, weighted=TRUE)
#' @rdname adjacencyMatrix-methods
#' @aliases adjacencyMatrix
#' @export
setMethod("adjacencyMatrix", signature(x="DNEA"),
          adjacencyMatrix.DNEA)

adjacencyMatrixReplace.DNEA <- function(x, weighted=FALSE, value){

  if(weighted){
    x@adjacency_matrix$weighted_adjacency <- value
  } else if(!weighted){
    x@adjacency_matrix$unweighted_adjacency <- value
  }
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("adjacencyMatrix", signature(x="DNEA"),
                 adjacencyMatrixReplace.DNEA)

adjacencyGraph.DNEA <- function(x, graph){

  x@consensus_clustering@adjacency_graphs[[graph]]
}
#' Retrieve the adjacency graph for the case, control,
#' or joint network
#'
#' The function  returns the adjacency graph made for
#' the case, control, or joint network constructed via
#' consensus clustering using \code{\link{clusterNet}}.
#'
#'
#' @param x A \code{\link{DNEA}}
#'
#' @param graph A character string indicating which of
#' the adjacency graphs to return. Values can be "joint_graph"
#' for the whole graph object, or one of the group values
#' returned by \code{\link{networkGroups}}.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{clusterNet}}
#' @returns An \code{\link{igraph}} graph object corresponding
#' to the specified adjacency graph.
#' @examples
#' #dnw is a \code{\link{DNEA}} object with the results
#' #generated for the example data accessed by running
#' #data(TEDDY) in the console. The workflow for this data
#' #can be found in the vignette accessed by running
#' #browseVignettes("DNEA") in the console.
#' data("dnw")
#'
#' adjacencyGraph(dnw, graph="DM:case")

#' @rdname adjacencyGraph-methods
#' @aliases adjacencyGraph
#' @export
setMethod("adjacencyGraph", signature(x="DNEA"),
          adjacencyGraph.DNEA)

adjacencyGraphReplace.DNEA <- function(x, graph, value){

  x@consensus_clustering@adjacency_graphs$graph <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("adjacencyGraph", signature(x="DNEA"),
                 adjacencyGraphReplace.DNEA)

adjacencyGraph.consensusClusteringResults <- function(x, graph){

  x@adjacency_graphs$graph
}
#' @rdname adjacencyGraph-methods
#' @aliases adjacencyGraph
setMethod("adjacencyGraph", signature(x="consensusClusteringResults"),
          adjacencyGraph.consensusClusteringResults)

adjacencyGraphReplace.consensusClusteringResults <- function(x, graph, value){

  x@adjacency_graphs$graph <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("adjacencyGraph", signature(x="consensusClusteringResults"),
                 adjacencyGraphReplace.consensusClusteringResults)

summary.consensusClusteringResults <- function(object){

  object@summary
}
#' Retrieve a summary of a consensusClusteringResults object
#'
#' The function takes as input a
#' \code{\link{consensusClusteringResults}} object and returns
#' a summary of the results of consensus clustering determined
#' via \code{\link{clusterNet}}.
#'
#' @param object A \code{\link{consensusClusteringResults}} object.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{clusterNet}}
#' @returns A data frame that corresponds to a summary of the results
#' of consensus clustering.
#' @rdname summary.consensuClusteringResults-methods
#' @aliases summary.consensusClusteringResults
#' @keywords internal
#' @noRd
setMethod("summary", signature(object="consensusClusteringResults"),
          summary.consensusClusteringResults)

consensus_clustering.DNEA <- function(x){

  x@consensus_clustering
}
#' @keywords internal
#' @noRd
setMethod("consensus_clustering", signature(x="DNEA"),
          consensus_clustering.DNEA)

CCsummary.DNEA <- function(x){

  summary(x@consensus_clustering)
}
#' Retrieves the summary results of consensus clustering
#'
#' The function takes as input a \code{\link{DNEA}} object and
#' returns a summary  of the results of consensus clustering
#' stored in the consensus_clustering slot as a
#' \code{\link{consensusClusteringResults}} object.
#'
#' @param x A \code{\link{DNEA}} object.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{clusterNet}}
#' @returns A data frame summary of the consensus clustering
#' results from DNEA.
#' @examples
#' #dnw is a \code{\link{DNEA}} object with the results
#' #generated for the example data accessed by running
#' #data(TEDDY) in the console. The workflow for this data
#' #can be found in the vignette accessed by running
#' #browseVignettes("DNEA") in the console.
#' data(dnw)
#'
#' CCsummary(dnw)
#' @rdname CCsummary-methods
#' @aliases CCsummary
#' @export
setMethod("CCsummary", signature(x="DNEA"),
          CCsummary.DNEA)

CCsummaryReplace.DNEA <- function(x, value){

  x@consensus_clustering@summary <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("CCsummary", signature(x="DNEA"),
                 CCsummaryReplace.DNEA)

subnetworkMembership.DNEA <- function(x){

  x@consensus_clustering@subnetwork_membership
}
#' Retrieve the subnetwork membership for each feature
#'
#' The function takes as input a \code{\link{DNEA}}
#' object and returns the results of consensus clustering
#' determined via \code{\link{clusterNet}}.
#'
#' @param x A \code{\link{DNEA}} object.
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{clusterNet}}
#' @returns A data frame that corresponds to the results of
#' consensus clustering.
#' @examples
#' #dnw is a \code{\link{DNEA}} object with the results
#' #generated for the example data accessed by running
#' #data(TEDDY) in the console. The workflow for this data
#' #can be found in the vignette accessed by running
#' #browseVignettes("DNEA") in the console.
#' data(dnw)
#'
#' subnetworkMembership(dnw)
#' @rdname subnetworkMembership-methods
#' @aliases subnetworkMembership
#' @export
setMethod("subnetworkMembership", signature(x="DNEA"),
          subnetworkMembership.DNEA)

subnetworkMembership.consensusClusteringResults <- function(x){

  x@subnetwork_membership
}
#' @rdname subnetworkMembership-methods
#' @aliases subnetworkMembership
setMethod("subnetworkMembership", signature(x="consensusClusteringResults"),
          subnetworkMembership.consensusClusteringResults)

subnetworkMembershipReplace.DNEA <- function(x, value){

  x@consensus_clustering@subnetwork_membership <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("subnetworkMembership", signature(x="DNEA"),
                 subnetworkMembershipReplace.DNEA)

netGSAresults.DNEA <- function(x){

  x@netGSA
}
#' Access the netGSA slot of a DNEA object
#'
#' The function takes as input a \code{\link{DNEA}}
#' object and returns a summary of the enrichment
#' analysis results stored in the netGSA slot.
#'
#' @param x A \code{\link{DNEA}} object.
#'
#' @author Christopher Patsalis
#' @seealso
#' \code{\link{runNetGSA}},
#' \code{\link[netgsa:NetGSA]{netgsa::NetGSA()}}
#' @returns A data frame of the results from \code{\link{runNetGSA}}.
#' @examples
#' #dnw is a \code{\link{DNEA}} object with the results
#' #generated for the example data accessed by running
#' #data(TEDDY) in the console. The workflow for this data
#' #can be found in the vignette accessed by running
#' #browseVignettes("DNEA") in the console.
#' data(dnw)
#'
#' netGSAresults(dnw)
#' @rdname netGSAresults-methods
#' @aliases netGSAresults
#' @export
setMethod("netGSAresults", signature(x="DNEA"),
          netGSAresults.DNEA)

netGSAresultsReplace.DNEA <- function(x, value){

  x@netGSA <- value
  validObject(x)
  x
}
#' @keywords internal
#' @noRd
setReplaceMethod("netGSAresults", signature(x="DNEA"),
                 netGSAresultsReplace.DNEA)
