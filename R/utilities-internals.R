#' splits the input data by the specified condition
#'
#' split_by_condition will separate the input expression matrix into a list
#' of matrices. Each matrix corresponds to the expression data for one
#' condition specified by condition_levels
#'
#' @param dat a matrix of expression data wherein the features are rows
#' and samples are columns
#' @param condition_levels A list or vector of the unique conditions
#' present in dat
#' @param condition_by_sample A list or vector of the condition value
#' corresponding to each sample.
#' @returns A list of expression matrices
#' @keywords internal
#' @noRd
split_by_condition <- function(dat, condition_levels, condition_by_sample){

  ##check input
  if(is.null(colnames(dat))) stop("dat must have column names!")
  if(is.null(rownames(dat))) stop("dat must have row names!")
  if(any(is.null(names(condition_by_sample)))){
    stop("each element in condition_by_sample must be named for its corresponding sample!")
  }
  if(!vector_compare(colnames(dat), names(condition_by_sample))) {
    stop("The provided conditions do not correspond to the samples in expression matrix!")
  }

  ##create key for separating the data by key and running diagnostic tests, feature DE calculations
  separated_conditions_data <- vector(mode='list', length=length(condition_levels))
  names(separated_conditions_data) <- condition_levels
  for(cond in condition_levels){
    separated_conditions_data[[cond]] <- dat[,condition_by_sample == cond]
  }
  return(separated_conditions_data)

}

#' @keywords internal
#' @noRd
vector_compare <- function(vec1, vec2){
  tryCatch({all(vec1 == vec2)},
           error= function(e){
             message(conditionMessage(e))
             FALSE},
           warning = function(w){
             stop(message("Vectors to compare are not the same length!"))
             message(conditionMessage(w))
             FALSE})
}
#' @keywords internal
#' @noRd
table_metadata_check <- function(dat,
                                 sample_names,
                                 feature_names){
  if(!vector_compare(colnames(dat), sample_names)){
    stop("sample names in input data do not match samples in object!")
  }
  if(!vector_compare(rownames(dat), feature_names)){
    message("Features in input data do not match features in object!")
  }

  if(!any(inherits(dat, "matrix") | inherits(dat, "data.frame"))) {
    stop('the input metadata should be of class "matrix" or "data.frame"!')
  }

  if(!all(apply(dat, 1, function(x) is.numeric(x)))){
    stop("input data should be a numeric matrix or dataframe!")
  }
  return(TRUE)
}

#' @keywords internal
#' @noRd
object_data_check <- function(dat,
                              sample_names,
                              feature_names,
                              group_labels){
  output <- list()
  for(i in seq(length(dat))){
    dat[[i]] <- dat[[i]][feature_names, sample_names]
    tmp <- split_by_condition(dat[[i]],
                              condition_levels = levels(group_labels),
                              condition_by_sample = group_labels)
    if(length(tmp) != 2) stop("DNEA requires two conditions!")
    if(!all(names(tmp) %in% levels(group_labels))){
      stop("input list of matrices not named after experimental groups!")
    }
    if(!all(colnames(tmp[[1]]) %in% names(group_labels)) &
       !all(colnames(tmp[[2]]) %in% names(group_labels))){
      stop("group labels names do not match sample names!")
    }
    for(x in seq(length(tmp))){
      tmp_samples <- names(group_labels)[group_labels == levels(group_labels)[x]]
      table_metadata_check(tmp[[x]],
                           sample_names = tmp_samples,
                           feature_names = feature_names)
    }
    tmp <- list(tmp)
    names(tmp) <- names(dat)[i]
    output <- append(output, tmp)
  }

  return(output)
}
