#' @include JSEM.R
#'
NULL
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' BICtune
#' @export
#' @import gdata
#' @import zoo
#' @import igraph
#' @import glasso
#' @import glmnet
#' @import corpcor
#' @import dplyr
#' @import foreach
#' @import parallel
#' @import doParallel


BICtune <- function(object, nCores = 4, main.seed = 101){
  dat <- lapply(object@Dataset_summary$scaled_separated_conditions,
                                                              function(d) t(scale(t(d))))
  ## Joint estimation
  n4cov <- max(sapply(dat, ncol))

  #Pre-define a range of lambda to select the tuning parameters using BIC.
  lambda.guo <- seq(0.01, 0.3, 0.02)*sqrt(log(object@Dataset_summary$num_features)/n4cov)

  #This needs to be more informative based on the data.
  trainX <- t(do.call(cbind, dat))
  trainY <- c(rep(1, ncol(object@Dataset_summary$scaled_separated_conditions[[1]])),
             rep(2, ncol(object@Dataset_summary$scaled_separated_conditions[[2]])))

  cat("BIC using Guo et al ... \n")

  #initialize parallel process
  cl <- parallel::makeCluster(nCores)
  registerDoParallel(cl)
  bic.guo <- vector(mode = "list", length = length(lambda.guo))
  parallel::clusterExport(cl = cl, varlist = c("lambda.guo", "trainX","trainY", "main.seed"), envir = environment())
  parallel::clusterEvalQ(cl = cl, c(library("MASS"),library("glasso"), set.seed(main.seed)))
  parallel::clusterExport(cl = cl, varlist = c("CGM_AHP_tune","CGM_AHP_train", "matTr"))



  bic.guo <- parallel::clusterMap(cl = cl, fun = 'CGM_AHP_tune', lambda = lambda.guo, MoreArgs = list(trainX = trainX,testX = trainX,model = trainY,BIC = TRUE,eta = 0.1))


  on.exit(stopCluster(cl))
  output <- vector(mode = 'list', length = 3)
  names(output) <- c('MinLambda','Lambda', 'bic.guo')
  output[["MinLambda"]] <- lambda.guo[which.min(sapply(bic.guo, function(a) a$BIC))]
  output[["Lambda"]] <- lambda.guo
  output[["bic.guo"]] <- bic.guo
  object@BIC <- output
  return(object)
}

#' StabilitySelection
#' @export
#' @import gdata
#' @import zoo
#' @import igraph
#' @import glasso
#' @import glmnet
#' @import corpcor
#' @import dplyr
#' @import foreach
#' @import parallel
#' @import doParallel
StabilitySelection <- function(object, UnevenGroups = FALSE, nreps = 50, nCores = 4, main.seed = 101){

  ##make sure all bic values are finite and remove those that are not
  print('start')
  tmp = sapply(object@BIC$bic.guo, function(a) a$BIC)
  print('tmp passed')
  if (max(is.infinite(tmp))==1){
    bic.guo <- object@BIC$bic.guo[is.finite(tmp)]
    lambda.guo <- object@BIC$Lambda[is.finite(tmp)]
    lastar.guo <- object@BIC$Lambda[which.min(sapply(object@BIC$bic.guo, function(a) a$BIC))]
    updated_bic <- vector(mode = 'list', length = 3)
    names(updated_bic) <-c('MinLambda','Lambda','bic.guo')
    updated_bic[['bic.guo']] <-bic.guo
    updated_bic[['MinLambda']] <- lastar.guo
    updated_bic[['Lambda']] <-lambda.guo
    object@BIC <- updated_bic
  } else {
    lastar.guo <- object@BIC$Lambda[which.min(sapply(object@BIC$bic.guo, function(a) a$BIC))]
  }

  print("move to list")
  ##stability selection, which requires lastar.guo from the previous step
  listX = lapply(as.list(object@Dataset_summary$scaled_separated_conditions),
                 function(d) scale(t(d)))
  print('test1')
  #stab.guo <- vector(mode = "list", length = length(nCores))
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  stab.guo <- as.list(1:nCores)

  parallel::clusterExport(cl = cl, varlist = c("CGM_AHP_tune","CGM_AHP_train","CGM_AHP_stabsel_subsample","CGM_AHP_stabsel", "matTr"))
  parallel::clusterExport(cl = cl, varlist = c("stab.guo", "lastar.guo","nreps","main.seed", "UnevenGroups"), envir = environment())
  parallel::clusterEvalQ(cl = cl, c(library("MASS"),library("glasso"), set.seed(main.seed)))
  #Use multiple nCores to run the function my.iter()
  # So in total we get nreps*nCores subsampling for stability selection.
  if (UnevenGroups){
    cat("Stability selection with additional subsampling ... \n")
    cat("Stability selection with Guo et al ... \n")
    stab.guo <- parallel::clusterMap(cl = cl,
                                     fun = 'CGM_AHP_stabsel_subsample',
                                     stab.guo = stab.guo,
                                     MoreArgs = list(X = listX, cnt = nreps, lastar = lastar.guo))
    object@Stable.Networks <- stab.guo
  } else {
    cat("Stability selection without additional subsampling ... \n")
    cat("Stability selection with Guo et al ... \n")
    stab.guo <- parallel::clusterMap(cl = cl,
                                     fun = 'CGM_AHP_stabsel',
                                     stab.guo = stab.guo,
                                     MoreArgs = list(X = listX, cnt = nreps, lastar = object@BIC$MinLambda))
  }
  on.exit(stopCluster(cl))
  object@Stable.Networks <- stab.guo
  return(object)
}












