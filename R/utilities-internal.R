#' splits the input data by the specified condition
#'
#' split_by_condition will separate the input expression matrix into a list of matrices. Each matrix
#' corresponds to the expression data for one condition specified by condition_levels
#'
#' @param dat a matrix of expression data wherein the features are rows and samples are columns
#' @param condition_levels A list or vector of the unique conditions present in dat
#' @param condition_by_sample A list or vector of the condition value corresponding to each
#'        sample.
#' @return A list of expression matrices
#' @keywords internal
#' @noRd
split_by_condition <- function(dat, condition_levels, condition_by_sample){

  #check input
  if(ncol(dat) != length(condition_by_sample)) stop("The provided conditions do not correspond to the samples in expression matrix!")

  #create key for separating the data by key and running diagnostic tests, feature DE calculations
  separated_conditions_data <- vector(mode = 'list', length = length(condition_levels))
  names(separated_conditions_data) <- condition_levels

  for(cond in condition_levels){
    separated_conditions_data[[cond]] <- dat[,condition_by_sample == cond]

  }
  return(separated_conditions_data)

}
