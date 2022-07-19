#'Set "DNEAobject" class
#'
#'@import methods
#'@noRd
#Internal Function to create DNEAobject
setClass(Class = "DNEAobject",
         slots = c(
           Project.Name = 'character',
           Assays = 'list',
           Metadata = 'list',
           Dataset_summary = 'list',
           Nodes = 'list',
           Edges = 'list',
           BIC = 'list',
           Adjacency.Matrices = 'list',
           Stable.Networks = 'list',
           Joint.Graph = 'list',
           NetGSA = 'list'
           )
)

setValidity("DNEAobject", function(object){
  if(length(object@Project.Name) > 1){
    "@Project.Name must be a character string"
  }
  for (i in length(object@Assays)){
    if(class(object@Assays[[i]]) != 'matrix'){
      "@Assays must be an expression matrix"
    }
    if(length(colnames(object@Assays[[i]])) != length(unique(colnames(object@Assays[[i]])))){
        "@Assays must be an expression matrix where each column is a unique feature."
      }
    if(!(is.numeric(object@Assays[[i]]))){
      "@Assays must be a matrix with numeric values."
    }
    # if(all(rownames(object@Assays[[i]]) != object@Metadata$Samples)){
    #   "Samples are out of order"
    # }
    # if(all(colnames(object@Assays[[i]]) != object@Metadata$clean_Feature_Names)){
    #   "Features are out of order"
    # }
    }
  if(class(object@Metadata$Samples) != 'character'){
    "@Metadata$Samples should be of class character"
  }
  if(class(object@Metadata$Features) != 'character'){
    " @Metadata$Features should be of class character"
  }
  if(class(object@Metadata$clean_Feature_Names) != 'character'){
    "@Metadata$clean_Feature_Names should be of class character"
  }
  if(class(object@Metadata$Condition) != 'factor' | length(levels(object@Metadata$Condition)) != 2){
    "@Metadata$Condition should be of class factor with 2 levels"
  }
})

















