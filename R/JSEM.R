#' maTr computes the trace of the matrix
#' @param z square numerical matrix
#' @return trace of the input matrix (scalar)
#' @noRd
matTr <- function(z){

  res <- sum(diag(z))

  return(res)
}

#' CGM_AHP_train will train a glasso model
#'
#' This function takes the data and corresponding condition, as well as a lambda value as input and
#' trains a glasso model, outputting the results. This code comes from Guo et al. (2011)
#'
#' @param trainX a matrix of expression data wherein the samples are rows and features are columns.
#' @param trainY a vector of corresponding conditions for samples in trainX
#' @param lambda_value the lambda hyperparameter to use in the glasso model
#' @param adaptive_weight default parameter used to calculate the penalty matrix
#' @param eta default parameter ??. Default is 0.01
#' @param limkappa default parameter that acts as the limit for the condition number of the sample cov.
#'        Default is 1e+6
#'
#' @returns The ouput of a glasso model
#'
#' @import zoo
#' @import glasso
#' @import glmnet
#' @import corpcor
#' @importFrom stats cov
#' @noRd
CGM_AHP_train <- function(
  trainX,
  trainY,
  lambda_value,
  adaptive_weight = array(1, c(length(unique(trainY)), ncol(trainX), ncol(trainX))),
  eta = 0.01,
  limkappa = 1e+6
  ){

  ## Set the general paramters
  K <- length(unique(trainY))
  p <- ncol(trainX)
  diff_value <- 1e+10
  count <- 0
  tol_value <- 0.01
  max_iter <- 30

  ## Set the optimization parameters
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
#'initialize static tuning variables to pass to CGM_AHP_tune
#'
#'@param trainX a matrix of expression data wherein the samples are rows and features are columns
#'@param model A list of corresponding conditions for the training_data
#'
#'@return a set of variables corresponding to the following variables regarding
#' the training_data: number of features (p), number of conditions (K), number of lambda
#' values being tested (N), an a vector containing the number of samples per condition (n),
#' an initialized vector for BIC scores (BIC_score), an initialized vector for likelihood (likelihood)
#'
#' @noRd
tune_init <- function(
    trainX, #training data
    model, #labels for models.
    lambda #a vector of supplied lambda values
){
  init_param <- vector(mode = 'list', length = 6)
  names(init_param) <- c('num_features','num_conditions','len_lambda','num_samples_by_cond', 'bic_score', 'likelihood')


  init_param[['num_features']] <- dim(trainX)[2]
  init_param[['num_conditions']] <- length(unique(model))
  init_param[['len_lambda']] <- length(lambda)
  init_param[['bic_score']] <- rep(0, init_param[['len_lambda']])
  init_param[['likelihood']] <- rep(0, init_param[['len_lambda']])

  num_samples_by_cond <- rep(0, init_param[['num_conditions']])
  for (k in 1:init_param[['num_conditions']]){
    num_samples_by_cond[k] <- length(which(model == k))
  }
  init_param[['num_samples_by_cond']] <- num_samples_by_cond
  return(init_param)
}

#' CGM_AHP_tune optimizes the tuning parameter for Guo et al. (2011) method in CGM_AHP_train
#'
#' This function takes the data and corresponding condition, as well as a list of lambda values
#' as input in order to optimize the lambda hyper-parameter.
#'
#'@param trainX a matrix of expression data wherein the samples are rows and features are columns.
#'@param testX a matrix of expression data wherein the samples are rows and features are columns.
#'       testX should be identical to trainX.
#'@param model a vector of corresponding conditions for samples in trainX
#'@param X a vector of lambda values tested to find the optimal hyper-parameter
#'@param BIC a boolean indicating whether or not to calculate the BIC score for each lambda.
#'       Default is FALSE.
#'@param eps A significance cut-off for thresholding network interactions.
#'       Default is 1e-06
#'@param eta default parameter ??. Default is 0.01
#'@param limkappa default parameter that acts as the limit for the condition number of the sample cov.
#'       Default is 1e+6
#'
#'@returns A list containing the BIC and liklihood score for each lambda parameter evaluated.
#'
#' @import zoo
#' @import glmnet
#' @import corpcor
#' @importFrom stats cov
#' @noRd
CGM_AHP_tune <- function(
    trainX,   # training data
    testX,    # test data
    model,    # labels for models.
    X,   # a vector of supplied lambda
    BIC=FALSE, # whether to compute the bic.score
    eps=1e-06,
    eta=0.01,
    limkappa = 1e+6
){

  num_features <- dim(trainX)[2]
  num_conditions <- length(unique(model))
  num_lambdas <- length(X)
  num_samples_by_cond <- rep(0, num_conditions)

  for (cond in 1:num_conditions){
    num_samples_by_cond[cond] <- length(which(model == cond))
  }

  bic_score <- 0
  likelihood <- 0


  Omega.hat <- CGM_AHP_train(trainX = trainX, trainY = model, lambda_value = X, eta = eta)$OMEGA

  for (k in 1:num_conditions){
    data <- testX[which(model == k), ]
    empcov <- cov(data)
    if (kappa(empcov) > limkappa){
      empcov = empcov + eta * diag(num_features)
    }
    likelihood = likelihood + matTr(empcov %*% Omega.hat[[k]]) - log(det(Omega.hat[[k]]))
  }

  if (BIC){
    for (k in 1:num_conditions){
      no.edge = sum(abs(Omega.hat[[k]])>eps) - num_features
      empcov <- cov(trainX[which(model==k),])
      bic_score = bic_score +  matTr(empcov %*% Omega.hat[[k]]) - log(det(Omega.hat[[k]])) + log(num_samples_by_cond[k]) * no.edge/(2*num_samples_by_cond[k])
    }
  }


  out <- list(BIC = bic_score, likelihood = likelihood)

  return(out)
}

#' stabsel_init initializes static tuning variables to pass to stability selection
#' (CGM_AHP_stabsel | CGM_AHP_stabsel_subsample)
#'
#' @param listX A list containing matrices of the expression data split by condition. The samples
#'        are in rows and the features are columns.
#' @param nreps The number of reps performed in stability selection
#'
#' @return a set of variables corresponding to the following:
#'         number of features (num_features), number of conditions (num_conditions),
#'         a vector containing sample numbers by condition (num_samples),
#'         number of samples in smallest condition (min_num_samples),
#'         a matrix of stability selection results (selection_matrix), and a matrix containing
#'         edge selection (edge_matrix)
#'         an initialized vector for BIC scores (BIC_score), an initialized vector for likelihood (likelihood)
#'
#' @noRd
stabsel_init <- function(
    listX, #The expression data split by condition
    nreps #number of reps performed for stability selection
){
  #initialize output list
  init_param <- vector(mode = 'list', length = 6)
  names(init_param) <- c('num_features',
                         'num_conditions',
                         'num_samples',
                         'min_num_samples',
                         'selection_matrix',
                         'edge_matrix')

  #set static variables for analysis
  if (class(listX) == 'list') {
    num_conditions = length(listX)
    num_features = ncol(listX[[1]])
    num_samples = unlist(lapply(listX, nrow))
    min_num_samples <- min(num_samples)
  } else{
    num_features <- ncol(listX)
    num_conditions <- 1
    num_samples <- nrow(listX)
    min_num_samples <- num_samples
  }

  selection_matrix = vector("list", num_conditions)
  names(selection_matrix) <- c('control','case')
  edge_matrix = vector("list", num_conditions)

  for (k in 1:num_conditions){
    selection_matrix[[k]] = matrix(0, num_features, num_features)
    edge_matrix[[k]] = matrix(0, num_features*(num_features-1)/2, nreps)
  }

  #set output results
  init_param[['num_features']] <- num_features
  init_param[['num_conditions']] <- num_conditions
  init_param[['num_samples']] <- num_samples
  init_param[['min_num_samples']] <- min_num_samples
  init_param[['selection_matrix']] <- selection_matrix
  init_param[['edge_matrix']] <- edge_matrix


  return(init_param)
}

#' CGM_AHP_stabsel performs stability selection WITHOUT additional subsampling
#'
#' This function will take as input the expression data and optimized lambda in order to randomly sample
#' the data to perform stability selection. The function is based on the method described by Guo et al.
#' (2011) and is designed to be run with some variation of the lapply family of functions.
#'
#' @param listX A list containing matrices of the expression data split by condition. The samples
#'        are in rows and the features are columns.
#' @param X A vector corresponding to the number of reps to be run in stability selection
#' @param lastar the optimized lambda parameter
#' @param eta default parameter ??. Default is 0.01
#' @param limkappa default parameter that acts as the limit for the condition number of the sample cov.
#'       Default is 1e+6
#' @param main.seed Sets the seed for random number generation. This ensures reproducibility.
#'
#' @return a precision matrix for the network corresponding to the input data created via that rep
#'
#' @import zoo
#' @import glmnet
#' @import corpcor
#' @noRd
CGM_AHP_stabsel <- function(listX,
                            X,
                            init_param,
                            lastar,
                            eta=0.01,
                            limkappa=1e+6,
                            main.seed = 101) {

  ########################################
  #**Set seed and initialize parameters**#
  ########################################

  #set the seed for reproducibility in random sampling
  set.seed(X*100+main.seed)

  X1 = vector("list", init_param[['num_conditions']])
  X2 = vector("list", init_param[['num_conditions']])

  model1 = NULL
  model2 = NULL

  selection_matrix <- init_param[['selection_matrix']]

  #############################################################
  #**sample expression data WITHOUT SUBSAMPLING for each rep**#
  #############################################################

  for (k in 1:init_param[['num_conditions']]){

    #create index to randomly sample half of the available samples
    ind1 = sample(x = seq(1, init_param[['num_samples']][[k]]), size = init_param[['num_samples']][[k]]/2, replace = FALSE)

    #create index to select the other half of samples not in ind.1
    ind2 = seq(1, init_param[['num_samples']][[k]])[match(seq(1, init_param[['num_samples']][[k]]), ind1, 0) == 0]

    #grab the first group of samples indicated by ind.1
    X1[[k]] = listX[[k]][ind1, ]

    #grab the remaining samples indicated by ind.2
    X2[[k]] = listX[[k]][ind2, ]

    model1 = c(model1, rep(k, length(ind1)))
    model2 = c(model2, rep(k, length(ind2)))
  }

  #########################################
  #**train model and concatenate results**#
  #########################################
  tmp1 = try(CGM_AHP_train(trainX=do.call(rbind, X1), trainY=model1, lambda_value=lastar, limkappa = limkappa, eta=eta))
  tmp2 = try(CGM_AHP_train(trainX=do.call(rbind, X2), trainY=model2, lambda_value=lastar, limkappa = limkappa, eta=eta))

  if (inherits(tmp1, "try-error") || inherits(tmp2, "try-error")){
    warning("There might be some error!")
    next;
  }

  for (k in 1:init_param[['num_conditions']]){
    mat1 = tmp1$OMEGA[[k]]
    mat1[which(abs(mat1)>1e-5)] = 1
    diag(mat1) = 0
    mat2 = tmp2$OMEGA[[k]]
    mat2[which(abs(mat2)>1e-5)] = 1
    diag(mat2) = 0
    selection_matrix[[k]] = selection_matrix[[k]] + mat1 + mat2
  }

  return(list(mat = selection_matrix, stab_sel_rep = X))
}

#' CGM_AHP_stabsel_subsample performs stability selection WITH additional subsampling
#'
#' This function will take as input the expression data and optimized lambda in order to randomly sample
#' the data to perform stability selection. The function is based on the method described by Guo et al.
#' (2011) and is designed to be run with some variation of the lapply family of functions. A key
#' difference between this function and the related CGM_AHP_stabsel is that random sampling is done
#' in a way that evens out sample distribution across the two conditions
#'
#' @param listX A list containing matrices of the expression data split by condition. The samples
#'        are in rows and the features are columns.
#' @param X A vector corresponding to the number of reps to be run in stability selection
#' @param lastar the optimized lambda parameter
#' @param eta default parameter ??. Default is 0.01
#' @param limkappa default parameter that acts as the limit for the condition number of the sample cov.
#'       Default is 1e+6
#' @param main.seed Sets the seed for random number generation. This ensures reproducibility.
#'
#' @return a precision matrix for the network corresponding to the input data created via that rep
#'
#' @import zoo
#' @import glmnet
#' @import corpcor
#' @noRd
CGM_AHP_stabsel_subsample <- function(listX,
                                      X,
                                      init_param,
                                      lastar,
                                      seed.base = 100,
                                      eta=0.01,
                                      limkappa=1e+6,
                                      main.seed = 101) {

  ########################################
  #**Set seed and initialize parameters**#
  ########################################

  #set the seed for reproducibility in random sampling
  set.seed(X*100+main.seed)

  #initialize selection matrix and edge matrix
  selection_matrix <- init_param[["selection_matrix"]]
  edge_matrix <- init_param[["edge_matrix"]]


  modelY = NULL
  subsampled_listX = vector("list", init_param[['num_conditions']])

  ##########################################################
  #**sample expression data WITH SUBSAMPLING for each rep**#
  ##########################################################

  ###subsampling
  subsampled_listX[[which.max(init_param[['num_samples']])]] = dplyr::sample_n(as.data.frame(listX[[which.max(init_param[['num_samples']])]]), 1.3*init_param[['min_num_samples']], replace = FALSE)

  temp90 = dplyr::sample_n(as.data.frame(listX[[which.min(init_param[['num_samples']])]]), 0.9*init_param[['min_num_samples']], replace = FALSE)
  temp10 = dplyr::sample_n(temp90, 0.1*init_param[['min_num_samples']], replace = FALSE)

  subsampled_listX[[which.min(init_param[['num_samples']])]] = rbind(temp90, temp10)

  #should we add a scaling step here?#

  subsampled_num_samples = lapply(subsampled_listX, nrow)

  for (k in 1:init_param[['num_conditions']]){
    modelY = c(modelY, rep(k, subsampled_num_samples[[k]]))
  }

  #########################################
  #**train model and concatenate results**#
  #########################################

  tmp_model = try(CGM_AHP_train(trainX=scale(do.call(rbind, subsampled_listX)), trainY=modelY, lambda_value=lastar, limkappa = limkappa, eta=eta))


  if (inherits(tmp_model, "try-error")){
    warning("There might be some error!")
    next;
  }

  for (k in 1:init_param[['num_conditions']]){
    tmp_mat = tmp_model$OMEGA[[k]]
    tmp_mat[which(abs(tmp_mat)>1e-5)] = 1
    diag(tmp_mat) = 0
    selection_matrix[[k]] = selection_matrix[[k]] + tmp_mat
    edge_matrix[[k]][,X] = t(tmp_mat)[lower.tri(tmp_mat,diag=F)]
  }

  return(list(mat = selection_matrix, edge_matrix = edge_matrix, stab_sel_rep = X))
}
#' adjDGlasso_minimal constructs the debiased glasso model
#'
#' This function takes the expression data and optimized lambda as input and fits the glasso model
#'
#' @param X A matrix of expression data wherein the samples are rows and features are columns.
#' @param weights A matrix of selection weights determined via stability selection to be integrated
#'        into the model. Default is 1.
#' @param theta_star The true precision matrix. Default is NULL
#' @param lambda The optimized lambda hyper-parameter. The default is NULL
#' @param quiet A boolean that indicates whether progress messages should be displayed. The default
#'        value is FALSE
#' @param zero.edge Indices of entries of inverse covariance to be constrained to zero (to be passed
#'        to glasso). The default is NULL
#'
#' @return An adjacency matrix for the data network estimated by the glasso model
#'
#' @import glasso
#' @importFrom stats cov2cor
#' @noRd
adjDGlasso_minimal <- function(
    X,
    weights=1,
    theta_star=NULL,
    lambda = NULL,
    quiet=FALSE,
    zero.edge=NULL
){
  n <- nrow(X)
  p <- ncol(X)
  X <- scale(X, center = T, scale = F)
  empcov <- (1/n) * (t(X) %*% X) #empirical cov
  if (is.null(lambda)){
    lambda <- sqrt(log(p)/n)
  }

  if (!is.null(zero.edge)){
    Theta.hat.from.Glasso <- glasso(s=empcov,
                                    rho=lambda*weights,
                                    penalize.diagonal=FALSE,
                                    zero=as.matrix(zero.edge))$wi
  } else {
    Theta.hat.from.Glasso <- glasso(s=empcov,
                                    rho=lambda*weights,
                                    penalize.diagonal=FALSE)$wi
  }
  Theta.hat.from.Glasso <- (Theta.hat.from.Glasso + t(Theta.hat.from.Glasso))/2

  if (!quiet){message("model estimated!\n", appendLF = TRUE)}

  coeff <- diag(1,p) - cov2cor(Theta.hat.from.Glasso)

  return(list(Theta.glasso=coeff))
}
