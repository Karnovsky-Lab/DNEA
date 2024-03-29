#' @include all-methods.R
NULL

#' Initialize DNEAobj object
#'
#' @description
#' This function takes as input a matrix of non-normalized, non-transformed
#' expression data and the case/control group labels in order to initiate a
#' DNEAobj object. Differential expression analysis using student's T-test and
#' Benjamini-Hochberg for multiple-testing corrections as well as diagnostic
#' testing are also performed on the data.
#'
#' ## IMPORTANT
#' Special attention should be given to the diagnostic criteria that is
#' output. The minimum eigen value and condition number are calculated for
#' the whole dataset as well as for each condition to determine mathematic
#' stability of the dataset and subsequent results from a GGM model. More
#' information about interpretation can be found in \strong{\emph{Details}}.
#'
#'
#' @param project_name A character string name for the experiment
#' @param expression_data A numeric matrix or data frame of un-transformed,
#' un-scaled expression data. The sample names should be row names and the
#' feature names should be column names
#' @param group_labels A factor vector of experimental group labels
#' named with the corresponding sample name
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{BICtune}},\code{\link{stabilitySelection}}
#'
#' @details
#' ## Diagnostics Motivation
#' Negative or zero eigenvalues in a dataset can represent
#' instability in that portion of the matrix, thereby invalidating
#' parametric statistical methods and creating unreliable results. In this
#' function, the minimum eigenvalue of the dataset is calculated by first
#' creating a pearson correlation matrix of the data. Instability may then
#' occur for a number of reasons, but one common cause is highly correlated
#' features (in the positive and negative direction). \cr
#'
#' Regularization often takes care of this problem by arbitrarily
#' selecting one of the variables in a highly correlated group and removing
#' the rest. We have developed DNEA to be very robust in situations where
#' \strong{\emph{p >> n}} by optimizing the model via several regularization
#' steps (\emph{please see} \code{\link{BICtune}} \emph{and}
#' \code{\link{stabilitySelection}}) that may handle such problems without
#' intervention, however, the user can also pre-emptively collapse
#' highly-correlated features into a single group via
#' \code{\link{aggregateFeatures}}.
#'
#' ## Benefits of Feature Aggregation
#' When your dataset contains highly correlated features, we recommend
#' aggregating features into related groups - such as highly-correlated
#' features of a given class of molecules (ie. many fatty acids,
#' carnitines, etc.) - because the user then has more control over which
#' variables are included in the model. Without collapsing, the model
#' regularization may result in one of the features within a class being
#' included and some or all of the remaining features being removed. By
#' collapsing first, you retain the signal from all of the features in the
#' collapsed group and also have information pertaining to which features
#' are highly correlated and will therefore have similar
#' feature-feature associations.
#'
#' @returns A \code{\link{DNEAobj}} object.
#'
#' @examples
#' #import example data
#' data(TEDDY)
#' data(T1Dmeta)
#' #create group labels
#' group_labels <- factor(T1Dmeta$group,
#'                        levels=c("DM:control", "DM:case"))
#' names(group_labels) <- rownames(T1Dmeta)
#'
#'
#' #initiate DNEAobj object
#' DNEA <- createDNEAobject(expression_data=TEDDY,
#'                          project_name="TEDDYmetabolomics",
#'                          group_labels=group_labels)
#'
#' @export
createDNEAobject <- function(project_name,
                             expression_data,
                             group_labels){


  ##restructure data to initiate object
  if(!missing(expression_data) & !missing(group_labels)){

    #check that order of group_labels matches data
    if(!(all(names(group_labels) == colnames(expression_data)))){
      stop("Group labels do not match sample order in expression data!")
    }

    ##turn condition into factor
    if(!is.factor(group_labels)){

      group_labels <- factor(group_labels)

      message("Condition for expression_data should be of class factor. ",
              "Converting now.\n",
              "Condition is now a factor with levels:\n1. ",
              levels(group_labels)[1],"\n2. ",
              levels(group_labels)[2])
    }

    #create data structures to initialize DNEAobject with un-scaled data
    restructured_data <- restructure_input_data(expression_data=expression_data,
                                                condition_values=group_labels)
  } else{
    #no data was provided - throw error
    stop('Expression data must be provided to create DNEAobject')
  }

  ##initiate DNEA object
  object <- new("DNEAobj",
                project_name=project_name,
                assays= restructured_data[[1]],
                metadata=list(samples=restructured_data[[2]]$samples,
                                features=restructured_data[[2]]$features,
                                network_group_IDs=structure(restructured_data[[2]]$samples$conditions,
                                                            names=restructured_data[[2]]$samples$samples),
                                network_groups=levels(restructured_data[[2]]$samples$conditions)),
                hyperparameter=list(BIC_scores=NULL, optimized_lambda=NULL, tested_lambda_values=NULL),
                adjacency_matrix=list(weighted_adjacency=NULL, unweighted_adjacency=NULL),
                stable_networks=list(selection_results=NULL, selection_probabilities=NULL))

  #check object
  validObject(object)

  #perform diagnostic testing on dataset
  datasetSummary(object) <- dataDiagnostics(mat=expressionData(x= object, assay="scaled_expression_data"),
                                       condition_values=networkGroups(object),
                                       conditions=networkGroupIDs(object))

  #print dataset summary
  message("\nDiagnostic criteria are as follows: ")
  show(object@dataset_summary)

  ##perform differential expression on the features
  nodeList(object) <- metabDE(mat=expressionData(x=object, assay="log_input_data"),
                              condition_values=networkGroups(object),
                              conditions=networkGroupIDs(object))

  #check for valid object once more
  validObject(object)

  return(object)
}
################################################################################
#' Restructure input data for initiation of DNEAobj object
#'
#' This function takes as input a matrix of expression data and the
#' experimental group labels in order to restructure the input as to prepare
#' it for initiation of a DNEAobj object.
#'
#' @param expression_data A matrix or data frame of expression data.
#' The sample names should be rownames and the feature names should be
#' column names. Column 1 should be a factor of the two conditions, followed by
#' the numeric expression data
#' @param condition_values A factor vector of experimental group labels named
#' with the corresponding sample name
#'
#' @returns A list containing two lists. The first list, named assays,
#' contains the uns-scaled data (if provided) in position 1 and the
#' scaled data in position 2. The second list, named metadata, contains
#' the metadata parsed from the input data.
#'
#' @author Christopher Patsalis
#'
#' @import methods
#' @importFrom janitor make_clean_names
#' @keywords internal
#' @noRd
restructure_input_data <- function(expression_data,
                                   condition_values){

  ##initialize output data structures
  #create metadata list and add names
  meta_key <- c("samples", "features")
  metadata <- vector(mode='list', length=length(meta_key))
  names(metadata) <- meta_key

  #create assays list of expression data
  assays_key <- c('input_data', 'log_input_data', 'scaled_expression_data')
  assays <- vector(mode='list', length=length(assays_key))
  names(assays) <- assays_key



  ##grab relevant metadata
  feature_names <- rownames(expression_data)
  clean_feature_names <- make_clean_names(feature_names)
  sample_names <- colnames(expression_data)

  ##convert expression data to matrix
  expression_data <- as.matrix(expression_data)

  ##clean column names to avoid R conflicts
  rownames(expression_data) <- clean_feature_names

  #add to assays list
  assays[['input_data']] <- expression_data
  assays[['log_input_data']] <- log(expression_data)

  ##scale the expression data for analysis
  #split by group
  scaled_expression_data <- lapply(split_by_condition(dat=expression_data,
                                                      condition_levels=levels(condition_values),
                                                      condition_by_sample=condition_values),
                                   function(x) t(scale(log(t(x)))))

  #order group_labels to make sure that still matches
  condition_values <- condition_values[colnames(expression_data)]

  message("Data has been normalized for further analysis.",
  " New data can be found in the scaled_expression_data assay!\n")

  ##concatenate output
  #add scaled data to assays list
  assays[['scaled_expression_data']] <- scaled_expression_data

  #Add features, samples, condition to metadata
  metadata[["samples"]] <- data.frame(samples=sample_names,
                                      conditions=condition_values,
                                      row.names=sample_names)

  metadata[["features"]] <- data.frame(feature_names=feature_names,
                                       clean_feature_names=clean_feature_names,
                                       row.names=feature_names)

  return(list(assays, metadata))
}

#' Calculate diagnostic criteria to determine stability of dataset
#'
#' This function takes as input a DNEAobj object and first splits the
#' scaled data by condition. The minimum eigen value and condition number
#' are then calculated for the whole dataset as well as for each condition to
#' determine mathematic stability of the dataset and subsequent results from
#' a GGM model. More information about interpretation can be found in
#' \strong{\emph{Details}}.
#'
#' @param mat A list of matrices corresponding to the normalized expression
#' data (one for the whole dataset and one for each condition)
#' @param condition_values A vector of condition names present
#' in the dataset
#' @param conditions A vector of condition labels corresponding to
#' the samples in mat
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{createDNEAobject}},\code{\link{aggregateFeatures}}
#'
#' @details
#' ## Diagnostics Motivation
#' Negative or zero eigenvalues in a dataset can represent instability
#' in that portion of the matrix, thereby invalidating parametric
#' statistical methods and creating unreliable results. In this function,
#' the minimum eigenvalue of the dataset is calculated by first creating
#' a pearson correlation matrix of the data. Instability may then occur
#' for a number of reasons, but one common cause is highly correlated
#' features (in the positive and negative direction). \cr
#'
#' Regularization often takes care of this problem by arbitrarily
#' selecting one of the variables in a highly correlated group and
#' removing the rest. We have developed DNEA to be very robust in
#' situations where \strong{\emph{p >>> n}} by optimizing the model via
#' several regularization steps (\emph{please see} \code{\link{BICtune}}
#' \emph{and} \code{\link{stabilitySelection}}) that may handle such
#' problems without intervention, however, the user can also pre-emptively
#' collapse highly-correlated features into a single group via
#' \code{\link{aggregateFeatures}}.
#'
#' ## Benefits of Feauture Collapsing
#' In scenarios like this we recommend collapsing highly correlated
#' features into a single group - particularly if the dataset contains
#' many highly-correlated features of a given class of molecules
#' (ie. many fatty acids, carnitines, etc.) - because the user then has
#' more control over which variables are included in the model. Without
#' collapsing, the model regularization may result in one of the features
#' within a class being included and some or all of the remaining features
#' being removed. By collapsing first, you retain the signal from all of
#' the features in the collapsed group and also have information pertaining
#' to which features are highly correlated and as a result track each other.
#'
#' @returns A DNEA object containing an initialized node_list. Differential
#' Expression is performed on the features if un-scaled data is provided.
#' The min eigen value and condition number is also printed for the whole
#' dataset as well as each condition.
#'
#' @importFrom stats t.test cor p.adjust
#' @keywords internal
#' @noRd
dataDiagnostics <- function(mat, condition_values, conditions) {

  ##set necessary parameters
  scaled_input_data <- do.call(cbind, mat)
  num_features <- nrow(scaled_input_data)
  num_samples <- ncol(scaled_input_data)

  ##min eigen value and condition number for each of the datasets
  #group1
  pearson_1 <- cor(t(mat[[condition_values[[1]]]]))
  cond_number_1 <- kappa(pearson_1, exact=TRUE)
  min_eigen_1 <- min(unlist(eigen(pearson_1, symmetric=TRUE, only.values=TRUE)))

  #group2
  pearson_2 <- cor(t(mat[[condition_values[[2]]]]))
  cond_number_2 <- kappa(pearson_2, exact=TRUE)
  min_eigen_2 <- min(unlist(eigen(pearson_2, symmetric=TRUE, only.values=TRUE)))

  #whole dataset
  pearson_3 <- cor(t(scaled_input_data))
  cond_number_3 <- kappa(pearson_3, exact=TRUE)
  min_eigen_3 <- min(unlist(eigen(pearson_3, symmetric=TRUE, only.values=TRUE)))

  ##create diagnostics table
  diagnostic_values <- data.frame(row.names=c("all_data",
                                              condition_values[[1]],
                                              condition_values[[2]]))
  diagnostic_values$min_eigen <- c(min_eigen_3,min_eigen_1, min_eigen_2)
  diagnostic_values$condition_num <- c(cond_number_3, cond_number_1, cond_number_2)

  if(any(diagnostic_values$min_eigen <= 1e-5)){

    warning("One or more condition(s) look unstable. ",
            "We recommend you aggregate features features before continuing!")
  } else{

    warning("Diagnostic tests complete. You can proceed with the analysis!")
  }

  return(new("DNEAinputSummary",
             num_samples=num_samples, num_features=num_features,
             diagnostic_values=diagnostic_values))

}

#' Perform differential expression analysis using students T-test
#'
#' Performs differential expression on features of the input data provided
#' using student's T-Test and multiple hypothesiscorrection
#' via Benjamini-Hochberg
#'
#' @param mat An m x n expression matrix with features as columns
#' @param condition_values a vector of the two possible condition values
#' @param conditions a vector containing the condition labels corresponding
#' to the samples in mat
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{createDNEAobject}}
#' @returns The differential expression results for the features in the
#' input matrix. The returned object is a data frame with the following
#' columns: \itemize{
#' \item \strong{clean feauture names -} The feature names
#' \item \strong{foldchange -} The foldchange in expression across
#' experimental conditions
#' \item \strong{fcdirection -} The direction of the foldchange,
#' ie. upregulated or downregulated
#' \item \strong{t-statistic -} The T-statistic from the student's
#' T-test used to perform Differential Expression
#' \item \strong{t-statistic direction -} The direction of the T-statistic
#' \item \strong{pvalue -} The associatated p-value for the differential
#' expression test
#' \item \strong{qvalue -} The p-value after adjusting for multiple-testing
#' using Benjamini-Hochberg
#' \item \strong{DEstatus -} Whether or not the features qvalue is less
#' than or equal to alpha=0.05}
#' @keywords internal
#' @noRd
metabDE <- function(mat,
                    condition_values,
                    conditions){

  ##set necessary parameters and output
  num_features <- nrow(mat)
  num_samples <- ncol(mat)
  feature_info <- data.frame(clean_feature_names=rownames(mat),
                             row.names=rownames(mat))

  ##split data by condition
  cond_data <- split_by_condition(dat=mat,
                                  condition_levels=condition_values,
                                  condition_by_sample=conditions)

  ##perform Differential expression
  feature_info$comparison <- rep(paste0(condition_values[2], ' / ', condition_values[1]), nrow(feature_info))
  feature_info$foldchange <- rowMeans(cond_data[[condition_values[2]]]) - rowMeans(cond_data[[condition_values[1]]])
  feature_info$fcdirection <- vapply(seq(1, num_features), function(i) ifelse(feature_info$foldchange[i] > 0, "Up", "Down"), FUN.VALUE=character(length(num_features)))

  metab_statistics <- vapply(rownames(feature_info), function(i) unlist(t.test(cond_data[[condition_values[2]]][i, ], cond_data[[condition_values[1]]][i, ], var.equal=FALSE)[c("statistic", "p.value")]), FUN.VALUE=numeric(2))
  feature_info$t_statistic <- metab_statistics["statistic.t",]
  feature_info$t_statistic_direction <- vapply(seq(nrow(feature_info)), function(i) ifelse(feature_info$t_statistic[i] > 0, "Up", "Down"), FUN.VALUE=character(1))
  feature_info$pvalue <- metab_statistics["p.value",]
  feature_info$qvalue <- p.adjust(feature_info$pvalue, "BH")
  feature_info$DEstatus <- vapply(seq(nrow(feature_info)), function(i) ifelse(abs(feature_info$qvalue[i]) >= 0.05, FALSE, TRUE), FUN.VALUE=logical(length(num_features)))

  return(feature_info)
}
