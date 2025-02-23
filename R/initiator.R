#' @include all-methods.R
NULL

#' Initialize a DNEA object
#'
#' @description
#' This function takes as input a matrix of non-normalized, non-transformed
#' expression data and the case/control group labels in order to initiate a
#' DNEA object. Differential expression analysis is performed using a
#' student's T-test and Benjamini-Hochberg for multiple-testing corrections.
#' Diagnostic testing is done on the input data by checking the minimum
#' eigen value and condition number of the expression data for each
#' experimental condition. To initialize a *DNEA* from a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment-class}},
#' or a mass_dataset-class from the massdataset package,
#' please see the \code{\link{sumExp2DNEA}} and
#' \code{\link{massDataset2DNEA}} documentation, respectively.
#'
#' ## IMPORTANT:
#' Special attention should be given to the diagnostic criteria that is
#' output. The minimum eigen value and condition number are calculated for
#' the whole data set as well as for each condition to determine mathematic
#' stability of the data set and subsequent results from a GGM model. More
#' information about interpretation can be found in the
#' \strong{\emph{Details}} section below.
#'
#'
#' @param project_name A character string name for the experiment.
#'
#' @param expression_data A numeric \emph{m x n} matrix or data frame
#' of un-transformed, un-scaled expression data. The sample names
#' should be column names and the feature names should be row names.
#'
#' @param scaled_expression_data A list of numeric \emph{m x n}
#' matrices or data frames of transformed and/or scaled expression
#' data. The sample names should be column names and the feature
#' names should be row names. Each set of expression data should
#' be aproximately normal.
#'
#' @param group_labels A factor vector of experimental group labels
#' named with the corresponding sample name.
#'
#' @param assay A character string indicating which assay to use for
#' diagnostics and differential expression analysis NOTE: The
#' function always defaults to using log transformed data for
#' differential expression analysis if provided.
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{BICtune}}, \code{\link{stabilitySelection}},
#' \code{\link{sumExp2DNEA}}, \code{\link{massDataset2DNEA}}
#'
#' @details
#' ## Diagnostics Motivation
#' Negative or zero eigenvalues in a data set can represent
#' instability in that portion of the matrix, thereby invalidating
#' parametric statistical methods and creating unreliable results. In this
#' function, the minimum eigenvalue of the data set is calculated by first
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
#' @returns A \code{\link{DNEA}} object.
#'
#' @examples
#' #import example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #create group labels
#' group_labels <- factor(T1Dmeta$group,
#'                        levels=c("DM:control", "DM:case"))
#' names(group_labels) <- rownames(T1Dmeta)
#'
#'
#' #initiate DNEA object
#' DNEA <- createDNEAobject(expression_data=TEDDY,
#'                          project_name="TEDDYmetabolomics",
#'                          group_labels=group_labels)
#' @export
createDNEAobject <- function(project_name,
                             expression_data,
                             scaled_expression_data,
                             group_labels,
                             assay){
  ##checks
  if(!is.character(project_name)) stop("project_name should be a string!")
  if(!missing(assay)){
    if(!is.character(assay)){
      stop("assay should be a string!")}}
  if(!missing(expression_data)){
    sample_names <- colnames(expression_data)
    feature_names <- rownames(expression_data)
  }
  if(!is.factor(group_labels)){
    group_labels <- factor(group_labels)
    message("Condition for expression_data should be of class factor. ",
            "Converting now.\n",
            "Condition is now a factor with levels:\n1. ",
            levels(group_labels)[1],"\n2. ",
            levels(group_labels)[2])
  }
  if(!missing(scaled_expression_data)){
    if(length(scaled_expression_data) != 1){
      stop("scaled_expression_data should be a named list containing",
           "a named list of expression data split by condition")
    }
    sample_names <- colnames(scaled_expression_data[[1]][[1]])
    feature_names <- rownames(scaled_expression_data[[1]][[1]])

    for(i in seq(length(scaled_expression_data))){
      names(scaled_expression_data[[i]]) <- NULL
      scaled_expression_data[[i]] <- do.call("cbind", scaled_expression_data[[i]])
      scaled_expression_data[[i]] <- scaled_expression_data[[i]][, names(group_labels)]
    }
  }
  ##restructure data
  if(!missing(expression_data)){
    restructured_data <- restructure_input_data(expression_data=expression_data,
                                                condition_values=group_labels)
    de_input <- restructured_data[["assays"]][["log_input_data"]]
    if(!missing(scaled_expression_data)){
      scaled_assays <- object_data_check(dat=scaled_expression_data,
                                         feature_names=feature_names,
                                         group_labels=group_labels)
      restructured_data[["assays"]] <- append(restructured_data[["assays"]],
                                              scaled_assays)
      if(missing(assay)){
        assay <- names(scaled_assays)[length(scaled_assays)]
      }
    }
    if(missing(assay)){assay <- "log-scaled_data"}
    ds_input <- restructured_data[["assays"]][[assay]]
  }else{
    scaled_assays <- object_data_check(dat=scaled_expression_data,
                                       feature_names=feature_names,
                                       group_labels=group_labels)
    restructured_data <- restructure_scaled_input_data(scaled_expression_data=scaled_assays,
                                                       condition_values=group_labels)
    warning("Raw peak intensity/concentration data is necessary for accurate analysis ",
            "in several steps of DNEA. If available, please input the raw ",
            "peak intensity/concentration data as well.")

    #data to be used in DE analysis later
    if(missing(assay)){
      assay <- names(scaled_assays)[length(scaled_assays)]
    }
    de_input <- do.call(cbind, restructured_data$assays[[assay]])
    ds_input <- restructured_data$assays[[assay]]
    warning("Scaling each condition prior to differential expression analysis may ",
            "lead to erroneous results. We strongly reccommend inputing ",
            "Raw peak intensity/concentration data!")
  }

  ##gather network information for use later
  network_group_IDs <- structure(restructured_data[["metadata"]]$samples$conditions,
                                 names=restructured_data[["metadata"]]$samples$samples)
  network_groups <- levels(restructured_data[["metadata"]]$samples$conditions)

  ##perform diagnostic testing on dataset
  message("Data diagnostics was performed on ", assay, " assay. ",
          "To check a different assay, please specify the assay parameter.")
  ds_test <- dataDiagnostics(mat=ds_input,
                             condition_values=network_groups,
                             conditions=network_group_IDs)

  ##perform differential expression on the features
  de_test <- metabDE(mat=de_input,
                     condition_values=network_groups,
                     conditions=network_group_IDs[colnames(de_input)])

  ##initiate DNEA object
  object <- new("DNEA",
                project_name=project_name,
                assays= restructured_data[["assays"]],
                metadata=list(samples=restructured_data[["metadata"]]$samples,
                              features=restructured_data[["metadata"]]$features,
                              network_group_IDs=network_group_IDs,
                              network_groups=network_groups),
                dataset_summary=ds_test, node_list=de_test,
                hyperparameter=list(BIC_scores=NULL, optimized_lambda=NULL,
                                    tested_lambda_values=NULL),
                adjacency_matrix=list(weighted_adjacency=NULL,
                                      unweighted_adjacency=NULL),
                stable_networks=list(selection_results=NULL,
                                     selection_probabilities=NULL))

  message("Diagnostic criteria are as follows: ")
  show(ds_test)
  return(object)
}

#' Initialize a DNEA from SummarizedExperiment
#'
#' @description
#' This function takes as input a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment-class}}
#'  object non-transformed in order to initiate a
#' \code{\link{DNEA}} object. Differential expression analysis
#' is performed using a student's T-test and Benjamini-Hochberg
#' for multiple-testing corrections. Diagnostic testing is done
#' on the input data by checking the minimum eigen value and
#' condition number of the expression data for each
#' experimental condition.
#'
#' ## IMPORTANT:
#' Special attention should be given to the diagnostic criteria that is
#' output. The minimum eigen value and condition number are calculated for
#' the whole data set as well as for each condition to determine mathematic
#' stability of the data set and subsequent results from a GGM model. More
#' information about interpretation can be found in the
#' \strong{\emph{Details}} section below.
#'
#'
#' @param object a SummarizedExperiment object
#'  object.
#'
#' @param scaled_expression_assay A character string corresponding to
#' the assay in the summarizedExperiment object to use for analysis.
#' Defaults to log-scaling the count data if not provided.
#'
#' @param group_label_col A character string corresponding to the
#' column in the sample metadata stored in the SummarizedExperiment object
#' to use as the group labels.
#' @inheritParams createDNEAobject
#'
#' @inherit createDNEAobject details
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{BICtune}}, \code{\link{stabilitySelection}},
#' \code{\link{createDNEAobject}}
#'
#' @returns A \code{\link{DNEA}} object.
#'
#' @examples
#' #load example data from airway package
#' library(airway)
#' data(airway)
#'
#' airway <- airway[1:50,]
#' airway <- airway[rowSums(SummarizedExperiment::assays(airway)$counts) > 5, ]
#' DNEA <- sumExp2DNEA(project_name = "airway",
#'                        object = airway,
#'                        group_label_col = "dex")
#'
#' @importFrom SummarizedExperiment assays colData rowData
#' @importFrom janitor make_clean_names
#' @export
sumExp2DNEA <- function(project_name,
                        object,
                        scaled_expression_assay,
                        group_label_col){
  ##checks
  if(!is.character(project_name)) stop("project_name should be a string!")
  if(!inherits(object, "SummarizedExperiment")){
    stop("object should inherit from class, SummarizedExperiment")
  }
  if(!(group_label_col %in% colnames(colData(object)))){
    stop("group_label_col should be a string corresponding to a",
         " column in colData(object)")
  }
  if(!missing(scaled_expression_assay)){
    if(!(scaled_expression_assay %in% names(assays(object)))){
      stop("scaled_expression_assay should be a string corresponding to",
           " the assay in assays(object) to use for analysis.")
    }
  }
  group_labels <- colData(object)[[group_label_col]]
  names(group_labels) <- rownames(colData(object))
  if(!is.factor(group_labels)){
    group_labels <- factor(group_labels)
    message("Condition for expression_data should be of class factor. ",
            "Converting now.\n",
            "Condition is now a factor with levels:\n1. ",
            levels(group_labels)[1],"\n2. ",
            levels(group_labels)[2])
  }
  expression_data <- SummarizedExperiment::assays(object)$counts + 1
  sample_names <- colnames(expression_data)
  feature_names <- rownames(expression_data)

  variable_info <- as.data.frame(rowData(object))
  sample_info <- as.data.frame(colData(object))
  tmp_dat <- as.list(SummarizedExperiment::assays(object))
  tmp_dat[["counts"]] <- tmp_dat[["counts"]] + 1
  tmp_dat <- object_data_check(dat = tmp_dat,
                               group_labels = group_labels,
                               sample_names = sample_names,
                               feature_names = feature_names)

  if(!missing(scaled_expression_assay)){
    scaled_data <- list()
    scaled_data[[scaled_expression_assay]] <- tmp_dat[[scaled_expression_assay]]
    output <- createDNEAobject(project_name = project_name,
                               expression_data = expression_data,
                               scaled_expression_data = scaled_data,
                               group_labels = group_labels)
  }else{
    scaled_expression_assay <- NULL
    output <- createDNEAobject(project_name = project_name,
                               expression_data = expression_data,
                               group_labels = group_labels)
  }

  #add sample metadata
  sample_info <- sample_info[rownames(metaData(output, type = "samples")),]

  output <- includeMetadata(object = output,
                            type = "samples",
                            metadata = sample_info)

  #add feature metadata
  variable_info <- variable_info[rownames(metaData(output, type = "features")),]
  output <- includeMetadata(object = output,
                            type = "features",
                            metadata = variable_info)

  #add additional assays
  new_assays <- names(tmp_dat)[-match(c("counts", scaled_expression_assay),
                                      names(tmp_dat))]
  for(i in new_assays){
    for(y in seq(length(tmp_dat))){
      rownames(tmp_dat[[i]][[y]]) <- make_clean_names(rownames(tmp_dat[[i]][[y]]))
    }
    output <- addExpressionData(object = output,
                                dat = tmp_dat[[i]],
                                assay_name = i)
  }

  return(output)
}
#' Initialize a DNEA object from a mass_dataset object
#'
#' @description
#' This function takes as input a
#' mass_dataset-class object from the massdataset package to
#' initiate a \code{\link{DNEA}} object. Differential
#' expression analysis is performed using a student's T-test
#' and Benjamini-Hochberg for multiple-testing corrections.
#' Diagnostic testing is done on the input data by checking
#' the minimum eigen value and condition number of the
#' expression data for each experimental condition.
#' \emph{\strong{NOTE: the massdataset package from the
#' tidymass software suite must be installed to use
#' this function. Please see
#' \url{https://massdataset.tidymass.org/} for more installation
#' instructions}}
#'
#' ## IMPORTANT:
#' Special attention should be given to the diagnostic criteria that is
#' output. The minimum eigen value and condition number are calculated for
#' the whole data set as well as for each condition to determine mathematic
#' stability of the data set and subsequent results from a GGM model. More
#' information about interpretation can be found in the
#' \strong{\emph{Details}} section below.
#'
#'
#' @param object a mass_dataset object.
#'
#' @param group_label_col A character string corresponding to the
#' column in the sample metadata stored in the mass_dataset object
#' to use as the group labels.
#'
#' @param scaled_input A TRUE/FALSE indicating whether the input data
#' is already normalized
#' @inheritParams createDNEAobject
#'
#' @inherit createDNEAobject details
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{BICtune}}, \code{\link{stabilitySelection}},
#' \code{\link{createDNEAobject}}
#'
#' @returns A \code{\link{DNEA}} object.
#'
#' @examples
#' #load data
#' data(TEDDY)
#' data(T1Dmeta)
#' data(metab_data)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#' T1Dmeta <- T1Dmeta[, c(6,7,7)]
#' colnames(T1Dmeta) <- c("sample_id", "group", "class")
#'
#' metab_data <- metab_data[rownames(TEDDY), ]
#'
#' sample_info_note = data.frame(name = c("sample_id", "group", "class"),
#'                               meaning = c("sample", "group", "class"))
#' variable_info_note = data.frame(name = c("variable_id", "mz", "rt"),
#'                                 meaning = c("variable_id", "mz", "rt"))
#' if (require(massdataset)) {
#' #create mass_dataset object from TEDDY
#' object <- massdataset::create_mass_dataset(expression_data = data.frame(TEDDY),
#'                                            sample_info = T1Dmeta,
#'                                            variable_info = metab_data,
#'                                            sample_info_note = sample_info_note,
#'                                            variable_info_note = variable_info_note)
#'
#' DNEA <- massDataset2DNEA(project_name = "mass_dataset",
#'                              object = object,
#'                              group_label_col = "group")
#' }
#'
#' @export
massDataset2DNEA <- function(project_name,
                             object,
                             group_label_col,
                             scaled_input=FALSE){
  ##checks
  if (!requireNamespace("massdataset", quietly=TRUE))
    stop("Could not load package massdataset. Is it installed?\n  ",
         "Note that massDataset2DNEA requires the tidymass package.\n  ",
         "Please install it with 'BiocManager::install(\"massdataset\")'.")
  if(!inherits(object, "mass_dataset")){
    stop("object should inherit of class mass_dataset from the",
         "massdataset package!")
  }
  if(!is.character(project_name)) stop("project_name should be a string!")
  expression_data <- as.data.frame(massdataset::extract_expression_data(object))
  sample_info <- data.frame(massdataset::extract_sample_info(object))
  rownames(sample_info) <- sample_info$sample_id
  if(!(group_label_col %in% colnames(sample_info))){
    stop("group_label_col should be a string corresponding to a",
         " column in massdataset::extract_sample_info(object)")
  }
  group_labels <- sample_info[[group_label_col]]
  names(group_labels) <- sample_info[["sample_id"]]
  group_labels <- group_labels[colnames(expression_data)]
  if(!is.factor(group_labels)){
    group_labels <- factor(group_labels)
    message("Condition for expression_data should be of class factor. ",
            "Converting now.\n",
            "Condition is now a factor with levels:\n1. ",
            levels(group_labels)[1],"\n2. ",
            levels(group_labels)[2])
  }

  variable_info <- data.frame(massdataset::extract_variable_info(object))
  var_cols <- colnames(variable_info) %in% c("variable_id",
                                             "mz", "rt",
                                             "fc", "p_value",
                                             "p_value_adjust")
  variable_info <- variable_info[, var_cols]
  sample_names <- colnames(expression_data)
  feature_names <- rownames(expression_data)
  ##restructure data
  tmp_dat <- list()
  tmp_dat[["mass_dataset"]] <- expression_data
  tmp_dat <- object_data_check(dat = tmp_dat,
                               group_labels = group_labels,
                               sample_names = sample_names,
                               feature_names = feature_names)
  if(scaled_input){
    output <- createDNEAobject(project_name = project_name,
                               scaled_expression_data = tmp_dat,
                               assay = names(tmp_dat),
                               group_labels = group_labels)
  }else{
    output <- createDNEAobject(project_name = project_name,
                               expression_data = expression_data,
                               group_labels = group_labels)
  }
  #add sample metadata
  sample_info <- sample_info[rownames(metaData(output, type = "samples")),]
  output <- includeMetadata(object = output,
                            type = "samples",
                            metadata = sample_info)

  #add feature metadata
  variable_info <- variable_info[rownames(metaData(output, type = "features")),]
  output <- includeMetadata(object = output,
                            type = "features",
                            metadata = variable_info)
  return(output)
}
################################################################################
#' Restructure input data for initiation of DNEA object
#'
#' This function takes as input a matrix of expression data and the
#' experimental group labels in order to restructure the input as to prepare
#' it for initiation of a \code{\link{DNEA}} object.
#'
#' @param expression_data A matrix or data frame of expression data.
#' The sample names should be column names and the feature names
#' should be row names.
#'
#' @param condition_values A factor vector of experimental group labels named
#' with the corresponding sample name
#'
#' @returns A list containing two lists. The first list, named "assays",
#' contains the un-scaled data (if provided) in position 1 and a list
#' of the scaled data for each condition in position 2.
#' The second list, named metadata, contains the metadata parsed
#' from the input data.
#'
#' @author Christopher Patsalis
#'
#' @import methods
#' @importFrom janitor make_clean_names
#' @keywords internal
#' @noRd
restructure_input_data <- function(expression_data,
                                   condition_values){
  if(!vector_compare(names(condition_values), colnames(expression_data))){
    stop("Group labels do not match sample order in expression data!")}

  ##initialize output data structures
  meta_key <- c("samples", "features")
  metadata <- vector(mode='list', length=length(meta_key))
  names(metadata) <- meta_key

  assays_key <- c('input_data', 'log_input_data', 'log-scaled_data')
  assays <- vector(mode='list', length=length(assays_key))
  names(assays) <- assays_key

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
  scaled_expression_data <- lapply(split_by_condition(dat=expression_data,
                                                      condition_levels=levels(condition_values),
                                                      condition_by_sample=condition_values),
                                   function(x) t(scale(log(t(x)))))

  ##order group_labels to make sure that still matches
  condition_values <- condition_values[colnames(expression_data)]

  message("Data has been normalized for further analysis.",
  " New data can be found in the log-scaled_data assay!")

  ##concatenate output
  assays[['log-scaled_data']] <- scaled_expression_data
  metadata[["samples"]] <- data.frame(samples=sample_names,
                                      conditions=condition_values,
                                      row.names=sample_names)
  metadata[["features"]] <- data.frame(feature_names=feature_names,
                                       clean_feature_names=clean_feature_names,
                                       row.names=feature_names)

  return(list(assays=assays, metadata=metadata))
}
################################################################################
#' Restructure scaled input data for initiation of DNEA object
#'
#' This function takes as input a list of matrices representing scaled
#' expression data and the experimental group labels in order to restructure
#' the input for initiation of a DNEA object.
#'
#' @param expression_data A list of matrices or data frames of expression
#' data. The sample names should be column names and the feature names
#' should be row names.
#'
#' @param condition_values A factor vector of experimental group labels named
#' with the corresponding sample name
#'
#' @returns A list containing two lists. The first list, named "assays",
#' contains the un-scaled data (if provided) in position 1 and a list
#' of the scaled data for each condition in position 2.
#' The second list, named metadata, contains the metadata parsed
#' from the input data.
#'
#' @author Christopher Patsalis
#'
#' @import methods
#' @importFrom janitor make_clean_names
#' @keywords internal
#' @noRd
restructure_scaled_input_data <- function(scaled_expression_data,
                                          condition_values){
  #compile all sample names
  sample_names <- c(names(condition_values)[condition_values == levels(condition_values)[[1]]],
                    names(condition_values)[condition_values == levels(condition_values)[[2]]])

  ##initialize output data structures
  meta_key <- c("samples", "features")
  metadata <- vector(mode='list', length=length(meta_key))
  names(metadata) <- meta_key

  assays_key <- names(scaled_expression_data)
  assays <- vector(mode='list', length=length(assays_key))
  names(assays) <- assays_key

  condition_values <- condition_values[sample_names]
  conditions_key <- levels(condition_values)

  feature_names <- rownames(scaled_expression_data[[assays_key[1]]][[conditions_key[1]]])
  clean_feature_names <- make_clean_names(feature_names)


  for(i in seq(length(scaled_expression_data))){
    ##convert expression data to matrix
    scaled_expression_data[[i]][[conditions_key[1]]] <- as.matrix(scaled_expression_data[[i]][[conditions_key[1]]])
    scaled_expression_data[[i]][[conditions_key[2]]] <- as.matrix(scaled_expression_data[[i]][[conditions_key[2]]])
    ##clean feature names to avoid R conflicts
    rownames(scaled_expression_data[[i]][[conditions_key[1]]]) <- clean_feature_names
    rownames(scaled_expression_data[[i]][[conditions_key[2]]]) <- clean_feature_names
    ##add to output
    assays[[assays_key[i]]] <- scaled_expression_data[[assays_key[i]]]
  }

  ##concatenate output
  metadata[["samples"]] <- data.frame(samples=sample_names,
                                      conditions=condition_values,
                                      row.names=sample_names)
  metadata[["features"]] <- data.frame(feature_names=feature_names,
                                       clean_feature_names=clean_feature_names,
                                       row.names=feature_names)

  return(list(assays=assays, metadata=metadata))
}



#' Calculate diagnostic criteria to determine stability of dataset
#'
#' This function takes as input a \code{\link{DNEA}} object and
#' calculates the minimum eigen value and condition number for the
#' whole data set as well as for each condition using the normalized
#' data. This is done to determine mathematic stability of the data
#' set and subsequent results from a GGM model. More information
#' about interpretation can be found in the
#' \strong{\emph{Details}} section.
#'
#' @param mat A list of matrices corresponding to the normalized
#' expression data accessed via \code{\link{expressionData}}.
#'
#' @param condition_values A vector of condition names present
#' in the data set.
#'
#' @param conditions A vector of condition labels corresponding to
#' the samples in mat.
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{createDNEAobject}},\code{\link{aggregateFeatures}},
#' \code{\link{expressionData}}
#'
#' @details
#' ## Diagnostics Motivation
#' Negative or zero eigenvalues in a data set can represent
#' instability in that portion of the matrix, thereby invalidating
#' parametric statistical methods and creating unreliable results. In this
#' function, the minimum eigenvalue of the data set is calculated by first
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
#' @returns A \code{\link{DNEA}} object containing an
#' initialized node_list. Differential expression is performed on the
#' features if un-scaled data is provided. The diagnostic criteria are
#' returned to the console.
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
  ##group1
  pearson_1 <- cor(t(mat[[condition_values[[1]]]]))
  cond_number_1 <- kappa(pearson_1, exact=TRUE)
  min_eigen_1 <- min(unlist(eigen(pearson_1, symmetric=TRUE, only.values=TRUE)))

  ##group2
  pearson_2 <- cor(t(mat[[condition_values[[2]]]]))
  cond_number_2 <- kappa(pearson_2, exact=TRUE)
  min_eigen_2 <- min(unlist(eigen(pearson_2, symmetric=TRUE, only.values=TRUE)))

  ##whole dataset
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

#' Perform differential expression analysis using a students T-test
#'
#' Performs differential expression on the features of the input data
#' using a student's T-Test and the Benjamini-Hochberg method for
#' multiple hypothesis correction
#'
#' @param mat An \emph{m x n} expression matrix.
#'
#' @param condition_values A vector of the two possible condition values.
#'
#' @param conditions A vector containing the condition labels corresponding
#' to the samples in mat.
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{createDNEAobject}}
#'
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
  feature_info <- data.frame(clean_feature_names=rownames(mat),
                             row.names=rownames(mat))
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
