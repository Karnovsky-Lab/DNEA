#library(dplyr)

#To compute the trace of a matrix
matTr <- function(z) sum(diag(z))


##--------------------------------------------\
## function: CGM_AHP_train.R
##--------------------------------------------\
## This is the code from Guo et al. (2011)
## trainX: data
## trainY: labels for categories (1, 2, 3,...)
## lambda_value: tuning parameter
## Make sure loading "glasso" package before using the code.
##--------------------------------------------\
CGM_AHP_train <- function(
  trainX, 
  trainY, 
  lambda_value,
  adaptive_weight = array(1, c(length(unique(trainY)), ncol(trainX), ncol(trainX))),
  eta = 0.01,
  limkappa = 1e+6,  #limit for condition number of the sample cov
  seed = 0,
  group = 0
){
  start.SS.time<-Sys.time()
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
  initialize.SS.time<-Sys.time()
  cat(paste('Group:', group, ', seed:', seed, ', CGM_AHP_train: maxiter is', max_iter, '. It has been', format(initialize.SS.time - start.SS.time, usetz = TRUE), 'from start of CGM_AHP_train. Parameters initialized and loop begins.'), sep = '\n', append = TRUE, file = Diagnostics)
  
  ## Start loop
  while ((count < max_iter) & (diff_value > tol_value)) {
    while.SS.time<-Sys.time()
    tmp <- apply(abs(OMEGA), c(2, 3), sum)
    tmp[abs(tmp) < 1e-10] <- 1e-10
    V <- 1/sqrt(tmp)
    
    for (k in seq(1, K)) {
      penalty_matrix <- lambda_value * adaptive_weight[k, , ] * V
      iter.SS.time<-Sys.time()
      cat(paste('Group:', group, ', seed:', seed, ', CGM_AHP_train: starting glasso with maxit 100 at', format(iter.SS.time-start.SS.time, usetz = TRUE), 'after CGM_AHP_train start.'),sep = '\n',append = TRUE, file = Diagnostics)
      obj_glasso <- glasso(S[k, , ], penalty_matrix, maxit = 100)
      for.end.SS.time<-Sys.time()
      cat(paste('Group:', group, ', seed:', seed, ', CGM_AHP_train: finished call to glasso. This call took', format(for.end.SS.time - iter.SS.time, usetz = TRUE), '. This call finished',
                format(for.end.SS.time - while.SS.time, usetz = TRUE), 'after the start of this iteration and', format(for.end.SS.time - start.SS.time, usetz = TRUE), 'after the start of CGM_AHP_train.'), sep = '\n',append = TRUE, file = Diagnostics)
      OMEGA_new[k, , ] <- (obj_glasso$wi + t(obj_glasso$wi))/2
      #OMEGA_new[k, , ] <- obj_glasso$wi
    }
    
    ## Check the convergence
    diff_value <- sum(abs(OMEGA_new - OMEGA))/sum(abs(OMEGA))
    count <- count + 1
    OMEGA <- OMEGA_new
    end.SS.time<-Sys.time()
    cat(paste('Group:', group, ', seed:', seed, ', CGM_AHP_train: count =', count, ', diff_value=', diff_value, '. This iter took', format(end.SS.time - while.SS.time, usetz = TRUE), 'to run and finished', format(end.SS.time - start.SS.time, usetz = TRUE), 'after CGM_AHP_train started.'), sep = '\n',append = TRUE, file = Diagnostics)
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
  
  func.SS.time<-Sys.time()
  cat(paste('CGM_AHP_train function took', format(func.SS.time - start.SS.time, usetz = TRUE), 'to complete.'), sep = '\n',append = TRUE, file = Diagnostics)
  rm(start.SS.time,initialize.SS.time,while.SS.time,iter.SS.time,for.end.SS.time,func.SS.time)
  return(output)
}



#to select the optimal tuning parameter for Guo et al.
##--------------------------------------------\
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
  start.SS.time<-Sys.time()
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
    iter.SS.time<-Sys.time()
    
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
    end.SS.time<-Sys.time()
    cat(paste("The", j, "-th step in tuning... This step took", format(end.SS.time - iter.SS.time, usetz = TRUE), 'and finished', format(end.SS.time - start.SS.time, usetz = TRUE), 'from the start of CGM_AHP_tune.'), sep = "\n", append = TRUE, file = Diagnostics)
  }
  
  out <- list(BIC = bic.score, likelihood = likelihood)
  func.SS.time<-Sys.time()
  cat(paste('CGM_AHP_tune function took', format(func.SS.time - start.SS.time, usetz = TRUE), 'to complete.'), sep = "\n", append = TRUE, file = Diagnostics)
  rm(start.SS.time, iter.SS.time, end.SS.time, func.SS.time)
  return(out)
}



##Stability selection for Guo's method conditional on lastar
##selected from cross validation
CGM_AHP_stabsel <- function(X, cnt, lastar, seed.base = 100, eta=0.01, limkappa=1e+6) {
  start.SS.time<-Sys.time()
  K = 1
  p = ncol(X)

  if (is.null(dim(X))) {
    K = length(X)
    p = ncol(X[[1]])
  }
  
  n = lapply(X, nrow)  
  
  X1 = vector("list", K)
  X2 = vector("list", K)
  sel.mat = vector("list", K)
  for (k in 1:K){
    sel.mat[[k]] = matrix(0, p, p)
  }
  count = 0
  initialize.SS.time<-Sys.time()
  cat(paste('CGM_AHP_stabsel function took', format(initialize.SS.time - start.SS.time, usetz = TRUE),'to initialize parameters.'), sep = "\n", append = TRUE, file = Diagnostics)
  for (i in 1:cnt) {
    iter.SS.time<-Sys.time()
    set.seed(i+seed.base)
    model.1 = NULL 
    model.2 = NULL 
    for (k in 1:K){
      ind.1 = sample(seq(1, n[[k]]), n[[k]]/2, F)
      ind.2 = seq(1, n[[k]])[match(seq(1, n[[k]]), ind.1, 0) == 0]
      X1[[k]] = X[[k]][ind.1, ]
      X2[[k]] = X[[k]][ind.2, ]
      model.1 = c(model.1, rep(k, length(ind.1)))
      model.2 = c(model.2, rep(k, length(ind.2)))
    }
    
    tmp.1 = try(CGM_AHP_train(trainX=do.call(rbind, X1), trainY=model.1, lambda_value=lastar, limkappa = limkappa, eta=eta, seed=i+seed.base, group=1))
    tmp.2 = try(CGM_AHP_train(trainX=do.call(rbind, X2), trainY=model.2, lambda_value=lastar, limkappa = limkappa, eta=eta, seed=i+seed.base, group=2))
    
    if (inherits(tmp.1, "try-error") || inherits(tmp.2, "try-error")){
      warning("There might be some error!")
      next;
    }
    
    for (k in 1:K){
      mat1 = tmp.1$OMEGA[[k]]
      mat1[which(abs(mat1)>1e-5)] = 1
      diag(mat1) = 0
      mat2 = tmp.2$OMEGA[[k]]
      mat2[which(abs(mat2)>1e-5)] = 1
      diag(mat2) = 0
      sel.mat[[k]] = sel.mat[[k]] + mat1 + mat2
    }
    
    count = count + 1
    end.SS.time<-Sys.time()
    cat(paste("CGM_AHP_stabsel: within-node rep count is:", i, '. This rep took', format(end.SS.time - iter.SS.time, usetz = TRUE),'. It has been', format(end.SS.time - start.SS.time, usetz = TRUE), 'since CGM_AHP_stabsel started.'), sep = '\n', append = TRUE, file = Diagnostics)
  }
  func.SS.time<-Sys.time()
  cat(paste('CGM_AHP_stabsel took', format(func.SS.time - start.SS.time, usetz = TRUE), 'to complete.'), sep = '\n', append = TRUE, file = Diagnostics)
  rm(start.SS.time, initialize.SS.time, iter.SS.time, end.SS.time, func.SS.time)
  return(list(mat = sel.mat, count = count))
}

##Stability selection for Guo's method conditional on lastar
##selected from cross validation
##subsampling
CGM_AHP_stabsel_subsample <- function(X, cnt, lastar, seed.base = 100, eta=0.01, limkappa=1e+6) {
  start.SS.time<-Sys.time()
  K = 1
  p = ncol(X)

  if (is.null(dim(X))) {
    K = length(X)
    p = ncol(X[[1]])
  }
  
  n = lapply(X, nrow)
  n_min = min(as.numeric(n))   
  
  sel.mat = vector("list", K)
  edge.mat = vector("list", K)
  
  for (k in 1:K){
    sel.mat[[k]] = matrix(0, p, p)
	edge.mat[[k]] = matrix(0, p*(p-1)/2, cnt)
  }
  count = 0
  initialize.SS.time<-Sys.time()
  cat(paste('CGM_AHP_stabsel_subsample parameters have been set and iterations begin. It has been',format(initialize.SS.time - start.SS.time, usetz = TRUE), 'since CGM_AHP_stabsel_subsample began.'), sep = '\n', append = TRUE, file = Diagnostics)
  for (i in 1:cnt) {
    iter.SS.time<- Sys.time()
    set.seed(i+seed.base)
	
	modelY = NULL
	new_X = vector("list", K)
	
	###subsampling
	new_X[[which.max(as.numeric(n))]] = dplyr::sample_n(as.data.frame(X[[which.max(as.numeric(n))]]), 1.3*n_min, replace = FALSE)
	
	temp90 = dplyr::sample_n(as.data.frame(X[[which.min(as.numeric(n))]]), 0.9*n_min, replace = FALSE)
	temp10 = dplyr::sample_n(temp90, 0.1*n_min, replace = FALSE)
	
	new_X[[which.min(as.numeric(n))]] = rbind.data.frame(temp90, temp10)
	
	#should we add a scaling step here?#
	
	new_n = lapply(new_X, nrow)
	
     for (k in 1:K){ 
	  modelY = c(modelY, rep(k, new_n[[k]]))
     }
    
	tmp.model = try(CGM_AHP_train(trainX=scale(do.call(rbind, new_X)), trainY=modelY, lambda_value=lastar, limkappa = limkappa, eta=eta))
		
	
    if (inherits(tmp.model, "try-error")){
      warning("There might be some error!")
      next;
    }
    
    for (k in 1:K){
      tmp.mat = tmp.model$OMEGA[[k]]
      tmp.mat[which(abs(tmp.mat)>1e-5)] = 1
      diag(tmp.mat) = 0
      sel.mat[[k]] = sel.mat[[k]] + tmp.mat
	  edge.mat[[k]][,i] = t(tmp.mat)[lower.tri(tmp.mat,diag=F)]
    }
    
    count = count + 1
    end.SS.time<-Sys.time()
    cat(paste("count is:", i, '. This loop took', format(end.SS.time - iter.SS.time, usetz = TRUE), 'and ended', format(end.SS.time - start.SS.time, usetz = TRUE), 'after start of CGM_AHP_stabsel_subsample.'), sep = '\n', append = TRUE, file = Diagnostics)
  }
  func.SS.time<-Sys.time()
  cat(paste("CGM_AHP_stabsel_subsample took", format(func.SS.time - start.SS.time, usetz = TRUE), 'to run.'), sep = '\n', append = TRUE, file = Diagnostics)
  rm(start.SS.time,initialize.SS.time,iter.SS.time,end.SS.time,func.SS.time)
  return(list(mat = sel.mat, count = count, edge.matrix = edge.mat))
}

##--------------------------------------------\
# Tuning parameter selection:
# Given a grid of lambda, compute the bic.score and
# likelihood based on a test dataset
##--------------------------------------------\
glasso_tune <- function(
  trainX, # training data n by p
  testX,  # test data n by p
  rho,    # a vector of tuning parameters
  weights, #the weight matrix for penalization
  eps=1e-06,
  eta=0.01, 
  limkappa=1e+6
){
  start.SS.time<-Sys.time()
  p = ncol(trainX)
  n = nrow(trainX)
  bic.score <- rep(0, length(rho))
  for (j in 1:length(rho)){
    iter.SS.time<-Sys.time()
    
    fit <- glasso(s = cov(trainX), rho = rho[j]*weights, penalize.diagonal = FALSE)
    sigInv = (fit$wi + t(fit$wi))/2
    sigInv[abs(sigInv)<eps] = 0
    
    empcov <- cov(testX) #empirical cov
    if (kappa(empcov) > limkappa){
      empcov = empcov + eta * diag(p)
    }
    bic.score[j] = matTr(empcov %*% fit$wi) - log(det(fit$wi)) +  + log(n) * (sum((abs(sigInv) > eps))-p)/(2*n)
    end.SS.time<-Sys.time()
    cat(paste("The ", j, "-th step in glasso tuning... This step took", format(end.SS.time - iter.SS.time, usetz = TRUE), 'and completed', format(end.SS.time - start.SS.time, usetz = TRUE), 'after the start of glasso_tune.'), sep = '\n', append = TRUE, file = Diagnostics)
  }
  
  out <- list(BIC = bic.score)
  func.SS.time<-Sys.time()
  cat(paste('glasso_tune took', format(func.SS.time - start.SS.time, usetz = TRUE), 'to run.'), sep = '\n', append = TRUE, file = Diagnostics)
  rm(start.SS.time, iter.SS.time, end.SS.time, func.SS.time)
  return(out)
}

require(glasso)
adjDGlasso <- function(
  X, #the n by p data matrix
  weights=1, #the weight for the penalty
  theta_star=NULL, #true precision matrix
  lambda = NULL, 
  FDR.type='BH', #FDR control procedure
  alpha = 0.1, #the significance level for FDR control
  quiet=TRUE,
  zero.edge=NULL # indices of entries of inverse covariance 
                   #to be constrained to be zero (to be passed to glasso)
){
  n <- nrow(X)
  p <- ncol(X)
  X <- scale(X, center = T, scale = F)
  empcov <- (1/n) * (t(X) %*% X) #empirical cov
  if (is.null(lambda)){
    lambda <- sqrt(log(p)/n)
  }
  
  if (!quiet){print('fit glasso')}
  
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
  
  if (!quiet){print('done')}
  
  if (!quiet){print('de-biasing glasso')}
  ## T.hat = Theta.hat - Theta.hat * (Sigma.hat - Theta.hat^{-1}) * Theta.hat
  ## Theta.hat and Sigma.hat are both symmetric.
  temp.mat <- empcov - chol2inv(chol(Theta.hat.from.Glasso))
  temp.vec <- as.vector(Theta.hat.from.Glasso %*% temp.mat %*% t(Theta.hat.from.Glasso))
  T.hat <- as.vector(Theta.hat.from.Glasso) - temp.vec
  T.hat <- matrix(T.hat,nrow = p)
  
  if (!quiet){print('done')}
  
  sigma.hat2 <- array(0,c(p,p))
  for (i in 1:p){
    for (j in 1:p){
      sigma.hat2[i,j] <- Theta.hat.from.Glasso[i,j]^2+Theta.hat.from.Glasso[i,i]*Theta.hat.from.Glasso[j,j]
    }
  }
  
  test.stat <- sqrt(n)*T.hat/sqrt(sigma.hat2)

  std.statistic <- NULL
  if (!is.null(theta_star)){
    std.statistic <- sqrt(n)*(T.hat - theta_star)/sqrt(sigma.hat2)
  }
  
  pvals <- 2*(pnorm(abs(test.stat), lower.tail=FALSE))
  pvals.vec <- lowerTriangle(pvals,diag=FALSE)
  adjpvals.vec <- p.adjust(pvals.vec, FDR.type)
  
  coeff <- diag(1,p) - cov2cor(T.hat)
  # coeff <- - pmax(pmin(T.hat, 1), -1)
  
  Qmat <- matrix(0, p, p)
  lowerTriangle(Qmat, diag=FALSE) <- adjpvals.vec
  Qmat <- Qmat + t(Qmat)
  diag(Qmat) <- rep(1, p)
  Qmat.fdr <- (Qmat <= alpha)
  
  # Qmat is the p by p matrix of qvalues;
  # Qmat.fdr is the thresholded matrix of qvalues based on alpha
  return(list(Theta=coeff, pvalue=pvals, qvalue = Qmat, qvalue.fdr=Qmat.fdr, statistic=std.statistic, theta.orig= Theta.hat.from.Glasso))
}

require(glasso)
adjDGlasso_minimal <- function(
  X, #the n by p data matrix
  weights=1, #the weight for the penalty
  theta_star=NULL, #true precision matrix
  lambda = NULL,
  quiet=FALSE,
  zero.edge=NULL # indices of entries of inverse covariance 
                   #to be constrained to be zero (to be passed to glasso)
){
  n <- nrow(X)
  p <- ncol(X)
  X <- scale(X, center = T, scale = F)
  empcov <- (1/n) * (t(X) %*% X) #empirical cov
  if (is.null(lambda)){
    lambda <- sqrt(log(p)/n)
  }
  
  if (!quiet){print('fit glasso')}
  
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
  
  if (!quiet){print('done')}

  coeff <- diag(1,p) - cov2cor(Theta.hat.from.Glasso)
  
  return(list(Theta.glasso=coeff))
}

##-----------------------------------\
##    StructDiff.full.R
##-----------------------------------\
### Purpose: to compare the estimated and the truth matrices. 
### Inputs: 
##    Ahat: a list of the estimated matrices. 
##    Amat: a list of the corresponding true matrices.
##  Outputs:
##   FPrate, FNrate, SHD, Floss, KL loss
##   The above measures are as defined in the paper, which are the average deviance
##   from all pairs of matrices.
##--------------------------------------------\
## The following is specifically for calculating the ROC curve
StructDiff <- function(Ahat, Amat, eps = 1e-08){
  K <- length(Amat)
  
  if (is.null(dim(Amat)) == F){
    # indicating there's only one matrix compared.
    K <- 1
    Ahat <- list(Ahat)
    Amat <- list(Amat)
  }
  
  p = dim(Amat[[1]])[1]
  TP <- rep(0, K)
  FP <- rep(0, K)
  TN <- rep(0, K)
  FN <- rep(0, K)
  SHD <- rep(0, K)
  FPrate <- rep(0, K)
  TPrate <- rep(0, K)
  FNrate <- rep(0, K)
  Pr <- rep(0, K)
  Re <- rep(0, K)
  F1 <- rep(0, K)
  Floss <- rep(0, K)
  for (k in 1:K){
    Floss[k] <- norm(Ahat[[k]] - Amat[[k]], "F")/norm(Amat[[k]], "F") 
    diag(Ahat[[k]]) = 0
    diag(Amat[[k]]) = 0
    
    TP[k] <- sum((abs(Ahat[[k]]) > eps) * (abs(Amat[[k]]) > eps))
    TN[k] <- sum((abs(Ahat[[k]]) <= eps) * (abs(Amat[[k]]) <= eps))
    FP[k] <- sum((abs(Ahat[[k]]) > eps) * (abs(Amat[[k]]) <= eps))
    FN[k] <- sum((abs(Ahat[[k]]) <= eps) * (abs(Amat[[k]]) > eps))  
    SHD[k] <- FP[k] + FN[k]
    
    P <- TP[k] + FN[k]
    N <- TN[k] + FP[k]
    TPrate[k] <- TP[k]/(P + eps)
    FPrate[k] <- FP[k]/(N + eps)
    FNrate[k] <- FN[k]/(P + eps)
    
    Re[k] <- TP[k]/(P + eps) ## Recall
    Pr[k] <- TP[k]/(TP[k] + FP[k] + eps) ## Precision
    
    F1[k] <- (2 * Pr[k] * Re[k])/(Pr[k] + Re[k] + eps)
    
  }
  
  dev = data.frame(TP = TP, FP = FP, TN = TN, FN = FN, 
                   TPrate = TPrate, FPrate = FPrate, FNrate = FNrate, 
                   SHD = SHD, Precision = Pr, Recall = Re, F1 = F1, FL = Floss)
  return(dev)
} 

