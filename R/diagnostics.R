#'Datadiagnostics
#'
#'@noRd
DATAdiagnostics <- function(object) {

  num_features <- length(featureNames(object))
  num_samples <- length(sampleNames(object))

  #Set up output structure
  feature_info <- data.frame('Features' = featureNames(object),
                             'clean_Feature_Names' = cleanFeatureNames(object))

  #create key for separating the data by key and running diagnostic tests, feature DE calculations
  separate_cond_key <- unlist(levels(condition(object)))


  #separate the dataset by condition
  if(!(is.null(Expression(object)))){

    separated_conditions_data <- vector(mode = 'list', length = length(levels(condition(object))))
    names(separated_conditions_data) <- separate_cond_key

    for(cond in separate_cond_key){
      separated_conditions_data[[cond]] <- t(Expression(object))[,(condition(object) == cond)]
    }
  }

  #separate the normalized data by condition
  separated_conditions_scaled <- vector(mode = 'list', length = length(levels(condition(object))))
  names(separated_conditions_scaled) <- separate_cond_key

  for(cond in separate_cond_key){
    separated_conditions_scaled[[cond]] <- t(NormalExpression(object))[,(condition(object) == cond)]
  }


  # if(!(is.null(NormalExpression(object)))){
  #
  #   separated_conditions_scaled <- vector(mode = 'list', length = length(levels(object@Metadata$Condition)))
  #   names(separated_conditions_scaled) <- separate_cond_key
  #   for(cond in separate_cond_key){
  #     separated_conditions_scaled[[cond]] <- t(object@Assays$NormalExpression)[,(object@Metadata$Condition == cond)]
  #   }
  #   if(is.null(object@Assays$Expression)){
  #     separated_conditions_data <-separated_conditions_scaled
  #   }
  # }

  #create list of conditions in data
  condition_levels <- vector(mode = 'list', length = length(levels(condition(object))))
  names(condition_levels) <- separate_cond_key
  condition_levels[[separate_cond_key[[1]]]]<-separate_cond_key[[1]]
  condition_levels[[separate_cond_key[[2]]]]<-separate_cond_key[[2]]

  # if(!(is.null(object@Assays$NormalExpression)) & is.null(object@Assays$Expression)){
  #   separated_conditions_scaled <- vector(mode = 'list', length = length(levels(object@Metadata$Condition)))
  #   names(separated_conditions_scaled) <- separate_cond_key
  #   for(cond in separate_cond_key){
  #     separated_conditions_scaled[[cond]] <- t(object@Assays$NormalExpression)[,(object@Metadata$Condition == cond)]
  #   }
  # } else{
  #   separated_conditions_scaled <- lapply(separated_conditions_data, function(d) t(scale(t(d))))
  # }
  #set data to use for feature information

  if(is.null(Expression(object))){
    cond_data <- separated_conditions_scaled
    separated_conditions_data <- NULL
    warning("No raw data was provided so Fold change calculation may be skewed")
  } else{
    cond_data <- separated_conditions_data
  }

  feature_info$foldchange <- rowMeans(cond_data[[2]]) - rowMeans(cond_data[[1]])
  feature_info$fcdirection <- sapply(1:num_features, function(i) ifelse(feature_info$foldchange[i] > 0, "Up", "Down"))
  feature_info$fc.notes <- paste0(levels(condition(object))[2], ' over ', levels(condition(object))[1])
  feature_info$statistic <- sapply(1:num_features, function(i) t.test(cond_data[[2]][i, ], cond_data[[1]][i, ], var.equal = FALSE)$statistic)
  feature_info$pvalue <- sapply(1:num_features, function(i) t.test(cond_data[[2]][i, ], cond_data[[1]][i, ], var.equal = FALSE)$p.value)
  feature_info$qvalue <- p.adjust(feature_info$pvalue, "BH")
  feature_info$DEstatus <- sapply(1:num_features, function(i) ifelse(abs(feature_info$qvalue[i]) >= 0.05, FALSE, TRUE))


  #Calculate important values
  pearson <- cor(scale(t(separated_conditions_scaled[[1]])))
  cond_number_1 <- kappa(pearson, exact=TRUE)
  min_eigen_1 <- min(unlist(eigen(pearson, symmetric=TRUE, only.values=TRUE)))
  pearson <- cor(scale(t(separated_conditions_scaled[[2]])))
  cond_number_2 <- kappa(pearson, exact=TRUE)
  min_eigen_2 <- min(unlist(eigen(pearson, symmetric=TRUE, only.values=TRUE)))
  pearson <- cor(scale(t(cbind(separated_conditions_scaled[[1]], separated_conditions_scaled[[2]]))))
  cond_number_3 <- kappa(pearson, exact=TRUE)
  min_eigen_3 <- min(unlist(eigen(pearson, symmetric=TRUE, only.values=TRUE)))

  Diagnostic_num <- data.frame(Dataset = c("Full_Dataset",unlist(levels(condition(object)))))
  Diagnostic_num$min_eigen <- list(min_eigen_3,min_eigen_1, min_eigen_2)
  Diagnostic_num$condition_num <- list(cond_number_3, cond_number_1, cond_number_2)

  #print values and suggestions
  cat('Diagnostic values are as follows:\n')
  print(as.matrix(Diagnostic_num))
  cat('\n')
  if(any(Diagnostic_num$min_eigen <=1e-5)){
    warning('One or more conditions look unstable. You should collapse features before continuing the analysis!\n')
  } else{
    cat('Diagnostic tests complete. You can proceed with the analysis!\n')
  }
  return(list(num_samples = num_samples, num_features = num_features,feature_summary = feature_info, Diagnostic_num = Diagnostic_num, condition_levels = condition_levels, separated_conditions = separated_conditions_data, scaled_separated_conditions = separated_conditions_scaled))
  #return(list(num_samples = num_samples, num_features = num_features,feature_summary = feature_info, Diagnostic_num = Diagnostic_num, condition_levels = condition_levels, separated_conditions = separated_conditions_data))

}




