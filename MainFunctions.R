createDNEAobject <- function(){
  dat <- read.csv(paste0(inFolder, filename,".csv"), header=TRUE, check.names = FALSE)
  colnames(dat)[2] <- "Sample_Group"

  dataset <- list()
  dataset$sample_info <- dat[,1:2]
  dataset$metab_info <- data.frame("Compound" =  colnames(dat)[-c(1,2)], stringsAsFactors = FALSE, check.names = FALSE)
  temp <- lapply(dat[,-c(1,2)], function(x) as.numeric(as.character(x)))
  names(temp) <- NULL
  dataset$dat <- do.call(rbind, temp)#p by n
  dataset$metab_info$ShortName <- dataset$metab_info$Compound
  print('Normalization Finished!')
  p <- nrow(dataset$metab_info)

  dat <- vector("list", ncond)
  dat[[1]] <- dataset$dat[,(dataset$sample_info$Sample_Group == group1)]
  dat[[2]] <- dataset$dat[,(dataset$sample_info$Sample_Group == group2)]

  dataset$metab_info$foldchange <- rowMeans(dat[[2]]) - rowMeans(dat[[1]])
  dataset$metab_info$fcdirection <- sapply(1:p, function(i) ifelse(dataset$metab_info$foldchange[i]>0, "Up", "Down"))
  dataset$metab_info$fc.notes <- "Group_2 over Group_1"

  dataset$metab_info$statistic <- sapply(1:p, function(i) t.test(dat[[2]][i,], dat[[1]][i,], var.equal=FALSE)$statistic)
  dataset$metab_info$pvalue <- sapply(1:p, function(i) t.test(dat[[2]][i,], dat[[1]][i,], var.equal=FALSE)$p.value)
  dataset$metab_info$qvalue <- p.adjust(dataset$metab_info$pvalue, "BH")
  dataset$metab_info$DEstatus <- sapply(1:p, function(i) ifelse(abs(dataset$metab_info$qvalue[i])>=0.05, FALSE, TRUE))

  save(dataset, file = paste0(OutFolder, filename,"_dataset_summary",".rda"))


  ## Joint estimation
  dat <- lapply(dat, function(d) t(scale(t(d))))
  n4cov <- max(sapply(dat, ncol))
  trainX <- t(do.call(cbind, dat))
  trainY <- c(rep(1, ncol(dat[[1]])), rep(2, ncol(dat[[2]])))

  lambda.guo = seq(0.01, 0.3, 0.02)*sqrt(log(p)/n4cov)

}


BICtune <- function(nCores = 4, nreps = 50){
  #Pre-define a range of lambda to select the tuning parameters using BIC.
  #This needs to be more informative based on the data.

  cat("BIC using Guo et al ... \n")
  cl <- makeCluster(nCores)
  registerDoParallel(cl)

  bic.guo = foreach(i = 1:length(lambda.guo),
                    .packages = c("MASS", "glasso")) %dopar%
    CGM_AHP_tune(trainX, testX=trainX, model=trainY, lambda=lambda.guo[i], BIC=TRUE, eta=0.1)

  stopCluster(cl)

  lastar.guo <- lambda.guo[which.min(sapply(bic.guo, function(a) a$BIC))]
  save(bic.guo, file=paste0(OutFolder, filename,"_BIC_tuning",".rda"))
}



StabilitySelection <- function(dat, UnevenGroups = FALSE, nreps = 50, nCores = 4, seed.base = 100){
  ##read in the necessary parameters from BIC tuning
  load(paste0(OutFolder,filename,"_BIC_tuning",".rda"))
  tmp = sapply(bic.guo, function(a) a$BIC)
  if (max(is.infinite(tmp))==1){
    bic.guo <- bic.guo[is.finite(tmp)]
    lambda.guo <- lambda.guo[is.finite(tmp)]
    lastar.guo <- lambda.guo[which.min(sapply(bic.guo, function(a) a$BIC))]
  } else {
    lastar.guo <- lambda.guo[which.min(sapply(bic.guo, function(a) a$BIC))]
  }
  ##stability selection, which requires lastar.guo from the previous step
  listX = lapply(dat, t)

  if (UnevenGroups){
    cat("Stability selection with additional subsampling ... \n")
    my.iter <- function(iter, seed.base){
      fit = CGM_AHP_stabsel_subsample(X=listX, cnt=nreps, lastar = lastar.guo, seed.base=seed.base)
      return(fit)
    }
  } else {
    cat("Stability selection without additional subsampling ... \n")
    my.iter <- function(iter, seed.base){
      fit = CGM_AHP_stabsel(X=listX, cnt=nreps, lastar = lastar.guo, seed.base=seed.base)
      return(fit)
    }
  }

  cat("Stability selection with Guo et al ... \n")

  #Use multiple nCores to run the function my.iter()
  # So in total we get nreps*nCores subsampling for stability selection.

  cl <- makeCluster(nCores)
  registerDoParallel(cl)

  stab_guo = foreach(i = 1:nCores,.packages = c("MASS", "glasso")) %dopar% my.iter(i,i*100+main.seed)

  stopCluster(cl)


  save(stab_guo, file=paste0(OutFolder,filename,"_stable_networks",".rda"))
}












