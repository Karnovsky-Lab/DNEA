#'
#'
#'
#'
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
#' @param scaled_expression_data A matrix or dataframe similar to expression_data but for scaled data
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
restructure_input_data <- function(expression_data, scaled_expression_data, control, case){

  ######################
  #**INITIALIZE LISTS**#
  #*####################

  #create metadata list and add names
  meta_key<-c("samples", "features", "clean_feature_names","condition_values")
  metadata<-vector(mode = 'list', length = length(meta_key))
  names(metadata) <- meta_key

  #create assays list of expression data
  assays_key <- c('expression_data','scaled_expression_data')
  assays <- vector(mode = 'list', length = length(assays_key))
  names(assays) <- assays_key


  ############################################################
  #**CHECK FOR PROPER CONDITIONS AND CONVERT DATA TO MATRIX**#
  ############################################################

  #un-scaled data
  if(missing(expression_data) == FALSE){

    #check proper data structure
    if(is.numeric(expression_data[,1])) stop('First column should be sample condition')

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

    #clean column names and convert data to matrix
    unscaled_feature_names <- colnames(expression_data[,-1])
    unscaled_sample_names <- rownames(expression_data)
    unscaled_condition_values <- expression_data[,1]
    colnames(expression_data) <- make_clean_names(colnames(expression_data))
    expression_data <- as.matrix(expression_data[,-1])

    #add to assays list
    assays[['expression_data']] <- expression_data
  }

  #scaled data
  if(missing(scaled_expression_data) == FALSE){

    #check proper data structure
    if(is.numeric(scaled_expression_data[,1])) stop('First column should be sample condition')

    #turn condition into a factor
    if(!(is.factor(scaled_expression_data[,1]))){
      scaled_expression_data[,1]<- factor(scaled_expression_data[,1],levels = c(control, case))
      cat(paste0('Condition for scaled_expression_data should be of class factor. Converting Now. \n',
                 'Condition is now a factor with levels:',
                 '\n', '1. ',
                 levels(scaled_expression_data[,1])[1],
                 '\n', '2. ',
                 levels(scaled_expression_data[,1])[2],'\n\n'))
    }

    #clean column names and convert data to matrix
    scaled_feature_names <- colnames(scaled_expression_data[,-1])
    scaled_sample_names <- rownames(scaled_expression_data)
    scaled_condition_values <- scaled_expression_data[,1]
    colnames(scaled_expression_data) <- make_clean_names(colnames(scaled_expression_data))
    scaled_expression_data <- as.matrix(scaled_expression_data[,-1])
  } else {

    #if no scaled data is provided, split by condition and scale the data
    scaled_expression_data <- lapply(split_by_condition(dat = expression_data,
                                                              condition_levels = levels(unscaled_condition_values),
                                                              condition_by_sample = unscaled_condition_values),
                                           function(x) scale(t(x)))

    #combine scaled data into one matrix
    scaled_expression_data <- rbind(scaled_expression_data[[1]],scaled_expression_data[[2]])

    #order samples to be same as un-scaled
    scaled_expression_data <- scaled_expression_data[rownames(expression_data),]

    message('Data has been normalized for further analysis. New data can be found in the scaled_expression_data assay!\n')

    #set the metadata equal to un-scaled counterpart
    scaled_feature_names <- unscaled_feature_names
    scaled_sample_names <- unscaled_sample_names
    scaled_condition_values <- unscaled_condition_values
  }

  ##############################################
  #**add data to initialized lists for output**#
  ##############################################

  #add scaled data to assays list
  assays[['scaled_expression_data']] <- scaled_expression_data

  #Add features, samples, condition to metadata
  metadata[["samples"]] <- data.frame(samples = scaled_sample_names,
                                      conditions = scaled_condition_values,
                                      row.names = scaled_sample_names)
  metadata[["features"]] <- data.frame(feature_names = scaled_feature_names,
                                       clean_feature_names = make_clean_names(scaled_feature_names),
                                       row.names = scaled_feature_names)

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
createDNEAobject <- function(project_name, expression_data, scaled_expression_data, control, case){

  ##check data input to make sure it matches
  if(missing(expression_data) == FALSE & missing(scaled_expression_data) == FALSE){
    if(all(colnames(expression_data[,-1]) != colnames(scaled_expression_data[,-1])) |
       all(rownames(expression_data[,-1]) != rownames(scaled_expression_data[,-1])) |
       all(expression_data[,1] != scaled_expression_data[,1])) stop("Data matrices must have identical features, samples, and condition values")
  }

  ##restructure data to initiate object
  if(missing(scaled_expression_data) == FALSE & missing(expression_data) == FALSE){

    #create data structures to initialize DNEAobject with scaled and un-scaled data
    restructured_data <- restructure_input_data(expression_data = expression_data,
                                                scaled_expression_data = scaled_expression_data,
                                                control = control,
                                                case = case)
  } else{

    if(missing(expression_data) == FALSE){

      #create data structures to initialize DNEAobject with un-scaled data
      restructured_data <- restructure_input_data(expression_data = expression_data,
                                                  control = control,
                                                  case = case)
    } else {

      if(missing(scaled_expression_data) == FALSE){

        #create data structures to initialize DNEAobject with scaled data
        restructured_data <- restructure_input_data(scaled_expression_data = scaled_expression_data,
                                                    control = control,
                                                    case = case)

      } else{

        #no data was provided - throw error
        stop('Expression data must be provided to create DNEAobject')
      }
    }

  }

  ##initiate DNEA object
  object <- new("DNEAobject",
                project_name = project_name,
                assays =  restructured_data[[1]],
                metadata = restructured_data[[2]],
                joint_graph = make_empty_graph(n = ncol(restructured_data[[1]][['scaled_expression_data']]),
                                               directed = TRUE),
                hyperparameter = list(BIC_scores = NULL, optimized_lambda = NULL, tested_lambda_values = NULL),
                adjacency_matrix = list(weighted_adjacency = NULL, unweighted_adjacency = NULL),
                stable_networks = list(selection_results = NULL, selection_probabilities = NULL))


   ##perform diagnostic testing on dataset
   diagnostic_values <- dataDiagnostics(mat = expressionData(object, type = "normalized"),
                                        condition_values = levels(conditions(object)),
                                        conditions = conditions(object))
   datasetSummary(object) <- diagnostic_values


   ##perform differential expression on the features
   DEresults <- metabDE(mat = expressionData(x = object, type = "input"),
                        condition_values = conditionLevels(object),
                        conditions = conditions(object))
   nodeList(object) <- DEresults

  return(object)
}
