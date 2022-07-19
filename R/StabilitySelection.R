StabilitySelection <- function(dat, UnevenGroups = FALSE, nreps = 50, nCores = 4, seed.base = 100){
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
