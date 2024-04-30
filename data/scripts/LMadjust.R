# LMadjust<- function (dataset){
#   dataset$Age<-unlist(lapply(dataset$Age, as.numeric))
#   adjusted_data<-dataset[,1:4]
#   for (i in 5:ncol(dataset)){
#     model_formula <-as.formula(paste(colnames(dataset)[i],'Age + Sex', sep = '~'))
#     lm.model<-lm(formula = model_formula, data = dataset)
#     new_data_points<-lm.model$residuals
#     new_data_points<-new_data_points + abs(min(new_data_points))
#     adjusted_data[[colnames(dataset)[i]]]<-with(adjusted_data, new_data_points)
#   }
#   adjusted_data <- log(adjusted_data[,5:ncol(adjusted_data)] + 1)
#   adjusted_data <- as.data.frame(scale(adjusted_data, center = TRUE, scale = TRUE))
#   output_data<-cbind.data.frame(dataset$sample,dataset$group,adjusted_data)
#   colnames(output_data)[1:2] <- c('sample', 'group')
#   return(output_data)
# }
#' LMadjust
#'
#' LMadjust takes in expression data and adjusts for the specified covariates
#'
#' LMadjust takes as input an expression table and a metadata table and uses linear regression to
#' adjust for the specified covariates. For each metabolite, a multiple linear regression is performed.
#' The residuals of the model are taken and the absolute value of the minimum residual +1 is added to all
#' of the residuals to create the new expression data.
#'
#' dataset - an expression table with rows as samples and columns as features. The row names should be
#' sample names and the column names metabolites
#'
#' covariates - A table of covariates for the data. The rows should be samples (with row names corresponding)
#' to the samples) and the rows are each covariate to be adjusted for.
#'
#' log_transform - a boolean whether or not to log transform
#' scale_data - a boolean whether or not to auto-scale data
LMadjust <- function(
    dataset,
    covariates,
    log_transform = FALSE,
    scale_data = FALSE){

  if(all(rownames(dataset) != rownames(covariates))) stop("covariates are in wrong order")

  cov_formula <- colnames(covariates)[1]
  for(cov in seq(2, length(covariates))){
    cov_formula<-paste(cov_formula,colnames(covariates[cov]), sep = ' + ')
  }

  dataset <- cbind(covariates,dataset)
  adjusted_data <- NULL
  for(i in seq((ncol(covariates) + 1), ncol(dataset))){
    model_formula <- as.formula(paste(colnames(dataset)[i],cov_formula, sep = '~'))
    lm.model<-lm(formula = model_formula, data = dataset)
    new_data_points<-lm.model$residuals
    new_data_points<-new_data_points + abs(min(new_data_points)) +1
    adjusted_data[[colnames(dataset)[i]]] <- new_data_points
    cat(paste('\r','column',i))
  }
  adjusted_data <- data.frame(adjusted_data)
  
  if(log_transform == TRUE){
    message('\ndata will be log transformed!')
    adjusted_data <- log(adjusted_data)
  } else{
    message('\ndata will be returned without log transforming!')
  }
  if(scale_data == TRUE){
    message("data will be scaled!")
    adjusted_data <- as.data.frame(scale(adjusted_data, center = TRUE, scale = TRUE))
  } else{
    message("data will be returned without scaling!")
  }

  return(adjusted_data)

}




