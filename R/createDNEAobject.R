#'createDNEAobject
#'@export
#'@import janitor
#'@import methods


restructure_input_data <- function(Expression = NULL, NormalExpression = NULL, control, case){

  #check to make sure data was input
  if(is.null(Expression) & is.null(NormalExpression)) stop('Expression data must be provided to create DNEAobject')

  #Check factor to make sure only 2 groups
  if(!(is.null(NormalExpression))){
    if(class(NormalExpression[,1]) == 'numeric') stop('First column should be sample condition')
    NormalExpression[,1]<- factor(NormalExpression[,1],levels = c(control, case))
    warning(paste0('Condition for NormalExpression should be of class factor. Converting Now. \n',
                   'Condition is now a factor with levels:','\n', '1. ', levels(NormalExpression[,1])[1], '\n', '2. ',levels(NormalExpression[,1])[2]))

  }
  if(!(is.null(Expression))){
    if(class(Expression[,1]) == 'numeric') stop('First column should be sample condition')
    if(class(Expression[,1]) != 'factor'){
      Expression[,1]<- factor(Expression[,1])
      warning(paste0('Condition for Expression should be of class factor. Converting Now. \n',
                     'Condition is now a factor with levels:','\n', '1. ', levels(Expression[,1])[1], '\n', '2. ',levels(Expression[,1])[2]))

    }
  }


  #set object to use for metadata
  if(!(is.null(Expression))){
    meta <- Expression
  } else{
    meta <- NormalExpression
  }


  #create metadata list and add names
  Meta_key<-c("Samples", "Features", "clean_Feature_Names","Condition")
  Metadata<-vector(mode = 'list', length = length(Meta_key))
  names(Metadata) <- Meta_key

  #Add features, samples, condition to metadata
  Metadata[["Samples"]] <- rownames(meta)
  Metadata[["Features"]] <- colnames(meta)[2:ncol(meta)]
  Metadata[["clean_Feature_Names"]] <- make_clean_names(colnames(meta)[2:ncol(meta)])
  Metadata[["Condition"]] <- meta[,1]

  #create Assays list of expression data
  Assays_key <- c('Expression','NormalExpression')
  Assays <- vector(mode = 'list', length = length(Assays_key))
  names(Assays) <- Assays_key

  #convert data to a matrix
  if(!(is.null(Expression))){
    colnames(Expression) <- make_clean_names(colnames(Expression))
    Expression <- as.matrix(Expression[,2:ncol(Expression)])

    # temp <- lapply(Expression[,-1], function(x) as.numeric(as.character(x)))
    # names(temp) <- NULL
    # Expression <- t(do.call(rbind, temp))#p by n
    # rownames(Expression)<-Metadata$Samples
    # colnames(Expression<-Metadata$clean_Feature_Names)

  }
  if(!(is.null(NormalExpression))){
    colnames(NormalExpression) <- make_clean_names(colnames(NormalExpression))
    NormalExpression <- as.matrix(NormalExpression[,2:ncol(NormalExpression)])
    # temp2 <- lapply(NormalExpression[,-1], function(x) as.numeric(as.character(x)))
    # names(temp2) <- NULL
    # NormalExpression <- t(do.call(rbind, temp2))
    # rownames(NormalExpression)<-Metadata$Samples
    # colnames(NormalExpression<-Metadata$clean_Feature_Names)
  }

  #Check to make sure data is the same
  if(!(is.null(Expression)) & !(is.null(NormalExpression))){
    #check to make sure the normalized and un-normalized input data have identical features and samples in the same order
    if(!(all(colnames(Expression) == colnames(NormalExpression))) | all(!(rownames(Expression) == rownames(NormalExpression)))) stop("Expression matrices must have identical features and samples")
  }

  #Normalize data if not provided
  if(is.null(NormalExpression)){
    NormalExpression <-scale(Expression)
    warning('Data has been normalized for further analysis. New data can be found in the NormalExpression Assay!')
  }

  #create matrices from data
  Assays[['Expression']] <- Expression
  Assays[['NormalExpression']] <- NormalExpression

  return(list(Assays, Metadata))

}
createDNEAobject <- function(Project.Name, Expression = NULL, NormalExpression = NULL, control, case){

  restructured_data <- restructure_input_data(Expression = Expression, NormalExpression = NormalExpression, control = control, case = case)

  #Initialize the DNEAobject
  object <- new("DNEAobject", Project.Name = Project.Name, Assays =  restructured_data[[1]], Metadata = restructured_data[[2]])
  #Perform diagnostics on the dataset
  object@Dataset_summary <- DATAdiagnostics(object)
  return(object)
}
