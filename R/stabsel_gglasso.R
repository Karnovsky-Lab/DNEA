# stability selection for group lasso

stabsel_gglasso <- function(X, Y, nSample, group_membership, method, seed.add){
	
	n <- nrow(X)
	half.n <- as.integer(n/2)
	p <- ncol(X)
	cluster.size <- table(group_membership)
	lambda.mat <- matrix(0, nrow = length(unique(group_membership)), 
							ncol = 2*nSample)
	
	for (i in seq(1,2*nSample,by=2)){
	  #cat("Count is: ", i, "\n")
		set.seed(i+seed.add)
		temp_seq <- sample(n)
		idx1 <- temp_seq[1:half.n]
		idx2 <- temp_seq[(half.n+1):n]
		gglasso1 <- gglasso(x = X[idx1,], y = Y[idx1], group = group_membership, loss = method)
		gglasso2 <- gglasso(x = X[idx2,], y = Y[idx2], group = group_membership, loss = method)
		
		lmbda1 <- numeric(length(cluster.size))
		lmbda2 <- numeric(length(cluster.size))
		lmbda1[1] <- round(gglasso1$lambda[max(which(gglasso1$beta[1,]==0))],2)
		lmbda2[1] <- round(gglasso2$lambda[max(which(gglasso2$beta[1,]==0))],2)
		
		for (a in 1:length(cluster.size)-1){
		
			temp1 <- round(gglasso1$lambda[max(which(gglasso1$beta[sum(cluster.size[1:a])+1,]==0))],2)
			temp2 <- round(gglasso2$lambda[max(which(gglasso2$beta[sum(cluster.size[1:a])+1,]==0))],2)
			lmbda1[a+1] <- temp1
			lmbda2[a+1] <- temp2
		}
		
		lambda.mat[,i] <- lmbda1
		lambda.mat[,i+1] <- lmbda2
		
	}
	
	return(lambda.mat)
}