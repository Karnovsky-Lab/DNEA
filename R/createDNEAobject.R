#'createDNEAobject
#'@export
#'@import janitor
#'@import methods
createDNEAobject <- function(Project.Name, Expression = NULL, NormalExpression = NULL, control, case, Metadata){

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

  }
  if(!(is.null(NormalExpression))){
    colnames(NormalExpression) <- make_clean_names(colnames(NormalExpression))
    NormalExpression <- as.matrix(NormalExpression[,2:ncol(NormalExpression)])
  }
  #Check to make sure data is the same
  if(!(is.null(Expression)) & !(is.null(NormalExpression))){
    #check to make sure the normalized and un-normalized input data have identical features and samples in the same order
    if(!(all(colnames(Expression) == colnames(NormalExpression))) | all(!(rownames(Expression) == rownames(NormalExpression)))) stop("Expression matrices must have identical features and samples")
    }

  #create matrices from data
  Assays[['Expression']] <- Expression
  Assays[['NormalExpression']] <- NormalExpression

  #Initialize the DNEAobject
  object <- new("DNEAobject", Project.Name = Project.Name, Assays =  Assays, Metadata = Metadata)

  #Perform diagnostics on the dataset
  object@Dataset_summary <- DATAdiagnostics(object)
  return(object)
}
