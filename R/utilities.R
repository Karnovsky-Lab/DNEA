#'
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

#' includeMetadata adds info to metadata
#'
#' This function will take additional metadata and add it to the corresponding dataframe in the metadata
#' slot.
#'
#' @param object A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#' @param type sample or feature metadata
#' @param metadata a dataframe containing metadata to add
#'
#' @return The same object as input with specified additions
#'
#' @export
includeMetadata <- function(object, type = c('sample', 'feature'), metadata){

  type = match.arg(type)
  if(type == 'sample'){
    if(all(sampleNames(object) == rownames(metadata))){
      for(i in 1:length(colnames(metadata))){


        object@metadata[["samples"]][[colnames(metadata)[i]]] <- metadata[, i]
      }
    } else{

      stop('new metadata order does not match sample order in DNEAobject')

    }
  } else{
    if(all(featureNames(object) == rownames(metadata)) |
       all(object@metadata[["features"]]$clean_feature_names == rownames(metadata))){
      for(i in 1:length(colnames(metadata))){

        new_metadata_colname <- colnames(metadata)[i]
        object@metadata[["features"]][[metadata[, i]]] <- metadata[, i]
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
#' @param object A pcorNetwork, DNEAobject, or DNEAobject_collapsed object
#'
#' @return The same object as input and saves the node and edge information as .csv files in
#'          the working directory
#' @importFrom utils write.csv
#' @export
getNetworkFiles <- function(object, file_path){

  if(missing(file_path)){
    file_path <- paste0(getwd(), "/")
  }
  #save node list
  write.csv(object@node_list, paste0(file_path, object@project_name,'_nodelist_',Sys.Date(),'.csv'),
            row.names = FALSE)

  #save edge list
  write.csv(object@edge_list, paste0(file_path, object@project_name,'_edgelist_',Sys.Date(),'.csv'),
            row.names = FALSE)

  return(object)
}

