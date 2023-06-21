#'
#'
#'
#'Structure input data for initialization of DNEAobject
#'
#'Manipulates input data and creates metadata for initialization of DNEAobject
#'
#' @param expression_data A matrix or dataframe of un-scaled expression data. The sample names should be rownames
#'        and the feature names should be column names. Column 1 should be a factor of the two conditions, followed by
#'        the numeric expression data
#' @param control A string corresponding to the name of condition 1 in column one of the data matrix
#' @param case A string corresponding to the name of condition 2 in column one of the data matrix
#'
#' @return A list containing two lists. The first list, named assays, contains the uns-scaled data (if provided)
#'         in position 1 and the scaled data in position 2. The second list, named metadata, contains the metadata parsed
#'         from the input data.
#'
#' @import methods
#' @importFrom janitor make_clean_names
#' @keywords internal
restructure_input_data <- function(expression_data, control, case){

  ######################
  #**INITIALIZE LISTS**#
  #*####################

  #create metadata list and add names
  meta_key<-c("samples", "features", "network_group_IDs", "network_groups")
  metadata<-vector(mode = 'list', length = length(meta_key))
  names(metadata) <- meta_key

  #create assays list of expression data
  assays_key <- c('expression_data','scaled_expression_data')
  assays <- vector(mode = 'list', length = length(assays_key))
  names(assays) <- assays_key


  ############################################################
  #**CHECK FOR PROPER CONDITIONS AND CONVERT DATA TO MATRIX**#
  ############################################################

  #check proper data structure
  if(!is.character(expression_data[,1])) stop('First column should be sample condition')

  #turn condition into factor
  if(!(is.factor(expression_data[,1]))){

    expression_data[,1]<- factor(expression_data[,1], levels = c(control, case))
    message(paste0('Condition for expression_data should be of class factor. Converting Now. \n',
                   'Condition is now a factor with levels:',
                   '\n', '1. ',
                   levels(expression_data[,1])[1],
                   '\n', '2. ',
                   levels(expression_data[,1])[2]))
  }

  #grab relevant metadata
  feature_names <- colnames(expression_data[,-1])
  clean_feature_names <- make_clean_names(feature_names)
  sample_names <- rownames(expression_data)
  condition_values <- expression_data[,1]

  #convert expression data to matrix
  expression_data <- as.matrix(expression_data[,-1])

  #clean column names to avoid R conflicts
  colnames(expression_data) <- clean_feature_names

  #add to assays list
  assays[['expression_data']] <- expression_data

  ##scale the expression data for analysis
  #split by group
  scaled_expression_data <- lapply(split_by_condition(dat = expression_data,
                                                      condition_levels = levels(condition_values),
                                                      condition_by_sample = condition_values),
                                   function(x) scale(t(x)))

  #combine scaled data into one matrix
  scaled_expression_data <- rbind(scaled_expression_data[[1]], scaled_expression_data[[2]])

  #order samples to be same as un-scaled
  scaled_expression_data <- scaled_expression_data[rownames(expression_data), ]

  message('Data has been normalized for further analysis. New data can be found in the scaled_expression_data assay!\n')

  ##############################################
  #**add data to initialized lists for output**#
  ##############################################

  #add scaled data to assays list
  assays[['scaled_expression_data']] <- scaled_expression_data

  #Add features, samples, condition to metadata
  metadata[["samples"]] <- data.frame(samples = sample_names,
                                      conditions = condition_values,
                                      row.names = sample_names)
  metadata[["features"]] <- data.frame(feature_names = feature_names,
                                       clean_feature_names = clean_feature_names,
                                       row.names = feature_names)

  return(list(assays, metadata))

}

#'Initializs the DNEAobject
#'
#'Takes the input data and restructures it to create the DNEAobject and expression_data is scaled if
#' scaled_exprssion_data is not provided. Diagnostics are performed on the scaled data and the minimum
#' eigen value and condition number are calculated for the whole dataset and each condition, specified by the
#' control and case inputs, separated.The output will inform the user if Feature Reduction should be performed prior
#' to continuing the analysis.
#'
#' @param project_name A string containing an identifying name for the analysis
#' @param expression_data A matrix or dataframe of un-scaled expression data. The sample names should be rownames
#'        and the feature names should be column names. Column 1 should be a factor of the two conditions, followed by
#'        the numeric expression data
#' @param scaled_expression_data A matrix or dataframe similar to expression_data but for scaled data
#' @param control A string corresponding to the name of condition 1 in column one of the data matrix
#' @param case A string corresponding to the name of condition 2 in column one of the data matrix
#'
#' @return DNEAobject containing the scaled and un-scaled (if provided) data in assays, as well as sample number,
#' feature number, sample names, feature names, and condition number in metadata
#'
#' @export
createDNEAobject <- function(project_name, expression_data, control, case){


  ##restructure data to initiate object
  if(!missing(expression_data)){

    #create data structures to initialize DNEAobject with un-scaled data
    restructured_data <- restructure_input_data(expression_data = expression_data,
                                                control = control,
                                                case = case)
  } else{

    #no data was provided - throw error
    stop('Expression data must be provided to create DNEAobject')
  }

  ##initiate DNEA object
  object <- new("DNEAresults",
                project_name = project_name,
                assays =  restructured_data[[1]],
                metadata = restructured_data[[2]],
                joint_graph = make_empty_graph(n = ncol(restructured_data[[1]][['scaled_expression_data']]),
                                               directed = TRUE),
                hyperparameter = list(BIC_scores = NULL, optimized_lambda = NULL, tested_lambda_values = NULL),
                adjacency_matrix = list(weighted_adjacency = NULL, unweighted_adjacency = NULL),
                stable_networks = list(selection_results = NULL, selection_probabilities = NULL))

  networkGroups(object) <- "conditions"

   ##perform diagnostic testing on dataset
   diagnostic_values <- dataDiagnostics(mat = expressionData(object, type = "normalized"),
                                        condition_values = networkGroups(object),
                                        conditions = networkGroupIDs(object))

   datasetSummary(object) <- new("DNEAinputSummary",
                                 num_samples = diagnostic_values[[1]],
                                 num_features = diagnostic_values[[2]],
                                 diagnostic_values = diagnostic_values[[4]])


   ##perform differential expression on the features
   DEresults <- metabDE(mat = expressionData(x = object, type = "input"),
                        condition_values = networkGroups(object),
                        conditions = networkGroupIDs(object))
   nodeList(object) <- DEresults

  return(object)
}
