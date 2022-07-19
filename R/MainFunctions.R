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


BICtune <- function(object, nCores = 4){
  ## Joint estimation
  # object@Dataset_summary$scaled_separated_conditions <- lapply(object@Dataset_summary$scaled_separated_conditions,
  #                                                              function(d) t(scale(t(d))))
  n4cov <- max(sapply(object@Dataset_summary$scaled_separated_conditions, ncol))

  #Pre-define a range of lambda to select the tuning parameters using BIC.
  lambda.guo <- seq(0.01, 0.3, 0.02)*sqrt(log(object@Dataset_summary$num_features)/n4cov)

  #This needs to be more informative based on the data.
  # trainX <- object@Assays$NormalExpression
  # testX <- object@Assays$NormalExpression
  # model <-unlist(lapply(object@Metadata$Condition, function(x) ifelse(x == object@Dataset_summary$condition_levels[[1]], 1,2)))
  trainX <- t(do.call(cbind, list(object@Dataset_summary$scaled_separated_conditions[[2]],object@Dataset_summary$scaled_separated_conditions[[1]])))
  testX <- t(do.call(cbind, list(object@Dataset_summary$scaled_separated_conditions[[2]],object@Dataset_summary$scaled_separated_conditions[[1]])))
  model <- c(rep(1, ncol(object@Dataset_summary$scaled_separated_conditions[[2]])),
             rep(2, ncol(object@Dataset_summary$scaled_separated_conditions[[1]])))
  BIC = TRUE
  eta = 0.1

  cat("BIC using Guo et al ... \n")

  #initialize parallel process
  cl <- parallel::makeCluster(nCores, type = "PSOCK")
  registerDoParallel(cl)
  bic.guo <- vector(mode = "list", length = length(lambda.guo))
  parallel::clusterEvalQ(cl = cl, c(library("MASS"),library("glasso"), set.seed(101)))
  parallel::clusterExport(cl = cl, varlist = c("CGM_AHP_tune","CGM_AHP_train", "matTr"))
  parallel::clusterExport(cl = cl, varlist = c("lambda.guo", "trainX", "testX","model", "BIC","eta"), envir = environment())


  bic.guo <- parallel::clusterMap(cl = cl, fun = 'CGM_AHP_tune', lambda = lambda.guo, MoreArgs = list(trainX,testX,model,BIC,eta))
  print('test2')

  on.exit(stopCluster(cl))
  output <- vector(mode = 'list', length = 3)
  names(output) <- c('MinLambda','Lambda', 'bic.guo')
  output[["MinLambda"]] <- lambda.guo[which.min(sapply(bic.guo, function(a) a$BIC))]
  output[["Lambda"]] <- lambda.guo
  output[["bic.guo"]] <- bic.guo
  object@BIC <- output
  return(object)
}

#' #To compute the trace of a matrix
#' #' maTr
#' #' @noRd
#' matTr <- function(z) sum(diag(z))


##--------------------------------------------\
## function: CGM_AHP_train.R
##--------------------------------------------\
## This is the code from Guo et al. (2011)
## trainX: data
## trainY: labels for categories (1, 2, 3,...)
## lambda_value: tuning parameter
## Make sure loading "glasso" package before using the code.
##--------------------------------------------\

#'CGM_AHP_train
#' @description CGM_AHP_train will
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
#' @noRd
CGM_AHP_train <- function(
    trainX,
    trainY,
    lambda_value,
    adaptive_weight = array(1, c(length(unique(trainY)), ncol(trainX), ncol(trainX))),
    eta = 0.01,
    limkappa = 1e+6  #limit for condition number of the sample cov
){
  ## Set the general paramters
  K <- length(unique(trainY))
  p <- ncol(trainX)
  diff_value <- 1e+10
  count <- 0
  tol_value <- 0.01
  max_iter <- 30

  ## Set the optimizaiton parameters
  OMEGA <- array(0, c(K, p, p))
  S <- array(0, c(K, p, p))
  OMEGA_new <- array(0, c(K, p, p))
  nk <- rep(0, K)

  ## Initialize Omega
  for (k in seq(1, K)) {
    idx <- which(trainY == k)
    S[k, , ] <- cov(trainX[idx, ])
    if (kappa(S[k, , ]) > limkappa) {
      S[k, , ] <- S[k, , ] + eta * diag(p)
    }
    tmp <- solve(S[k, , ])
    OMEGA[k, , ] <- tmp
    nk[k] <- length(idx)
  }

  ## Start loop
  while ((count < max_iter) & (diff_value > tol_value)) {
    tmp <- apply(abs(OMEGA), c(2, 3), sum)
    tmp[abs(tmp) < 1e-10] <- 1e-10
    V <- 1/sqrt(tmp)

    for (k in seq(1, K)) {
      penalty_matrix <- lambda_value * adaptive_weight[k, , ] * V
      obj_glasso <- glasso(S[k, , ], penalty_matrix, maxit = 100)
      OMEGA_new[k, , ] <- (obj_glasso$wi + t(obj_glasso$wi))/2
      #OMEGA_new[k, , ] <- obj_glasso$wi
    }

    ## Check the convergence
    diff_value <- sum(abs(OMEGA_new - OMEGA))/sum(abs(OMEGA))
    count <- count + 1
    OMEGA <- OMEGA_new
    #cat(count, ', diff_value=', diff_value, '\n')
  }

  ## Filter the noise
  list.Omega <- NULL
  Theta <- vector("list", K)
  for (k in seq(1, K)) {
    ome <- OMEGA[k, , ]
    ww <- diag(ome)
    ww[abs(ww) < 1e-10] <- 1e-10
    ww <- diag(1/sqrt(ww))
    tmp <- ww %*% ome %*% ww
    ome[abs(tmp) < 1e-05] <- 0
    OMEGA[k, , ] <- ome
    list.Omega[[k]] <- ome
    Theta[[k]] <- diag(diag(list.Omega[[k]])^(-0.5)) %*% list.Omega[[k]] %*% diag(diag(list.Omega[[k]])^(-0.5))
  }

  output <- list()
  output$OMEGA <- list.Omega
  output$S <- S
  output$Theta <- Theta
  output$lambda <- lambda_value

  return(output)
}


#to select the optimal tuning parameter for Guo et al.
##--------------------------------------------\
#' CGM_AHP_tune
#' @description CGM_AHP_tune will
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
#' @noRd
CGM_AHP_tune <- function(
    trainX,   # training data
    testX,    # test data
    model,    # labels for models.
    lambda,   # a vector of supplied lambda
    BIC=FALSE, # whether to compute the bic.score
    eps=1e-06,
    eta=0.01,
    limkappa = 1e+6
){
  p = dim(trainX)[2]
  K = length(unique(model))
  N <- length(lambda)
  n <- rep(0, K)
  for (k in 1:K){
    n[k] <- length(which(model == k))
  }

  bic.score <- rep(0, N)
  likelihood <- rep(0, N)
  for (j in 1:N){
    cat("The ", j, "-th step in tuning... \n")

    Omega.hat <- CGM_AHP_train(trainX = trainX, trainY = model, lambda_value = lambda[j], eta = eta)$OMEGA

    for (k in 1:K){
      data <- testX[which(model == k), ]
      empcov <- cov(data)
      if (kappa(empcov) > limkappa){
        empcov = empcov + eta * diag(p)
      }
      likelihood[j] = likelihood[j] + matTr(empcov %*% Omega.hat[[k]]) - log(det(Omega.hat[[k]]))
    }

    if (BIC){
      for (k in 1:K){
        no.edge = sum(abs(Omega.hat[[k]])>eps) - p
        empcov <- cov(trainX[which(model==k),])
        bic.score[j] = bic.score[j] +  matTr(empcov %*% Omega.hat[[k]]) - log(det(Omega.hat[[k]])) + log(n[k]) * no.edge/(2*n[k])
      }
    }
  }

  out <- list(BIC = bic.score, likelihood = likelihood)

  return(out)
}
#
# StabilitySelection <- function(dat, UnevenGroups = FALSE, nreps = 50, nCores = 4, seed.base = 100){
#   ##read in the necessary parameters from BIC tuning
#   load(paste0(OutFolder,filename,"_BIC_tuning",".rda"))
#   tmp = sapply(bic.guo, function(a) a$BIC)
#   if (max(is.infinite(tmp))==1){
#     bic.guo <- bic.guo[is.finite(tmp)]
#     lambda.guo <- lambda.guo[is.finite(tmp)]
#     lastar.guo <- lambda.guo[which.min(sapply(bic.guo, function(a) a$BIC))]
#   } else {
#     lastar.guo <- lambda.guo[which.min(sapply(bic.guo, function(a) a$BIC))]
#   }
#   ##stability selection, which requires lastar.guo from the previous step
#   listX = lapply(dat, t)
#
#   if (UnevenGroups){
#     cat("Stability selection with additional subsampling ... \n")
#     my.iter <- function(iter, seed.base){
#       fit = CGM_AHP_stabsel_subsample(X=listX, cnt=nreps, lastar = lastar.guo, seed.base=seed.base)
#       return(fit)
#     }
#   } else {
#     cat("Stability selection without additional subsampling ... \n")
#     my.iter <- function(iter, seed.base){
#       fit = CGM_AHP_stabsel(X=listX, cnt=nreps, lastar = lastar.guo, seed.base=seed.base)
#       return(fit)
#     }
#   }
#
#   cat("Stability selection with Guo et al ... \n")
#
#   #Use multiple nCores to run the function my.iter()
#   # So in total we get nreps*nCores subsampling for stability selection.
#
#   cl <- makeCluster(nCores)
#   registerDoParallel(cl)
#
#   stab_guo = foreach(i = 1:nCores,.packages = c("MASS", "glasso")) %dopar% my.iter(i,i*100+main.seed)
#
#   stopCluster(cl)
#
#
#   save(stab_guo, file=paste0(OutFolder,filename,"_stable_networks",".rda"))
# }












