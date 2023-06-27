my_test <- function(X){
  x <- Matrix(0, X, X)
  return(x)
}

my_test2 <- function(X){
  pkgs <- search()
  pkgs <- pkgs[grep("package:",pkgs)]
  y <- unlist(sapply(pkgs,lsf.str))
  return(list(pkgs = pkgs, func = y))
  }
cl<-parallel::makePSOCKcluster(4)
parallel::clusterExport(cl = cl, varlist = c("my_test"))
parallel::clusterEvalQ(cl = cl, c(library('Matrix')))
test.out <- parallel::parLapply(cl, fun = "my_test", X = 8:12)
parallel::stopCluster(cl)

my_test(7)

get('hey')
