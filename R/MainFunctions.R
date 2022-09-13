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
#' @import parallel


BICtune <- function(object, nCores = 4, main.seed = 101){

  # dat <- lapply(object@Dataset_summary$scaled_separated_conditions,
  #                                                             function(d) t(scale(t(d))))
  dat<- lapply(object@Dataset_summary$scaled_separated_conditions, function(d) t(scale(t(d))))

  ## Joint estimation
  n4cov <- max(sapply(dat, ncol))

  #Pre-define a range of lambda to select the tuning parameters using BIC.
  lambda.guo <- seq(0.01, 0.3, 0.02)*sqrt(log(numFeatures(object))/n4cov)

  #This needs to be more informative based on the data.
  trainX <- t(do.call(cbind, dat))
  trainY <- c(rep(1, ncol(object@Dataset_summary$scaled_separated_conditions[[1]])),
             rep(2, ncol(object@Dataset_summary$scaled_separated_conditions[[2]])))

  cat("BIC using Guo et al ... \n")

  #initialize parallel process
  cl <- parallel::makeCluster(nCores)
  #registerDoParallel(cl)
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
#' @import parallel
#' @import progress
StabilitySelection <- function(object, UnevenGroups = FALSE, nreps = 50, nCores = 4, main.seed = 101){

  ##make sure all bic values are finite and remove those that are not
  tmp = sapply(object@BIC$bic.guo, function(a) a$BIC)
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

  ##stability selection, which requires lastar.guo from the previous step
  listX = lapply(as.list(object@Dataset_summary$scaled_separated_conditions),
                 function(d) scale(t(d)))

  #stab.guo <- vector(mode = "list", length = length(nCores))

  #create independent processes to run reps in parallel
  cl <- makeCluster(nCores, outfile = "")
  #registerDoParallel(cl)

  stab.guo <- as.list(1:nCores)
  #create progress bar
  #p <- progressr::progressor(along = stab.guo)

  parallel::clusterExport(cl = cl, varlist = c("CGM_AHP_tune","CGM_AHP_train","CGM_AHP_stabsel_subsample","CGM_AHP_stabsel", "matTr"))
  parallel::clusterExport(cl = cl, varlist = c("stab.guo", "listX", "lastar.guo","nreps","main.seed", "UnevenGroups"), envir = environment())
  parallel::clusterEvalQ(cl = cl, c(library("MASS"),library("glasso"),library("progress"), set.seed(main.seed)))
  #Use multiple nCores to run the function my.iter()
  # So in total we get nreps*nCores subsampling for stability selection.
  if (UnevenGroups){
    cat("Stability selection with additional subsampling ... \n")
    cat("Stability selection with Guo et al ... \n")
    # stab.guo <- parallel::clusterMap(cl = cl,
    #                                  fun = 'CGM_AHP_stabsel_subsample',
    #                                  stab.guo = stab.guo,
    #                                  MoreArgs = list(X = listX, cnt = nreps, lastar = lastar.guo))
    stab.guo <- parallel::parLapply(cl = cl,
                                    fun ="CGM_AHP_stabsel_subsample",
                                    X = stab.guo,
                                    listX = listX,
                                    cnt = nreps,
                                    lastar = lastar.guo)
    object@Stable.Networks <- stab.guo
  } else {
    cat("Stability selection without additional subsampling ... \n")
    cat("Stability selection with Guo et al ... \n")
    # stab.guo <- parallel::clusterMap(cl = cl,
    #                                  fun = 'CGM_AHP_stabsel',
    #                                  stab.guo = stab.guo,
    #                                  MoreArgs = list(X = listX, cnt = nreps, lastar = object@BIC$MinLambda))
    stab.guo <- parallel::parLapply(cl = cl,
                                   fun ="CGM_AHP_stabsel",
                                   X = stab.guo,
                                   listX = listX,
                                   cnt = nreps,
                                   lastar = lastar.guo)
  }
  # if (UnevenGroups){
  #   cat("Stability selection with additional subsampling ... \n")
  #   cat("Stability selection with Guo et al ... \n")
  #   stab.guo <- pbmapply(cl = cl,
  #                       FUN = CGM_AHP_stabsel_subsample,
  #                       X = 1:nCores,
  #                       MoreArgs = list(listX = listX, cnt = nreps, lastar = lastar.guo))
  #   object@Stable.Networks <- stab.guo
  # } else {
  #   cat("Stability selection without additional subsampling ... \n")
  #   cat("Stability selection with Guo et al ... \n")
  #   stab.guo <- pbMapply(cl = cl,
  #                       FUN = CGM_AHP_stabsel,
  #                       X = 1:nCores,
  #                       MoreArgs = list(listX = listX, cnt = nreps, lastar = lastar.guo))
  # }
  on.exit(stopCluster(cl))
  object@Stable.Networks <- stab.guo
  return(object)
}
#' GetNetworks
#' @export
#' @import gdata
#' @import zoo
#' @import igraph
#' @import glasso
#' @import glmnet
#' @import corpcor
#' @import dplyr
GetNeworks <- function(object, nCores = 4, nreps = 2, UnevenGroups = FALSE, eps = 1e-06){

    ## Retrieve stable networks, which requires the stability selection results from previous step
    sel_mat <- vector("list", length(object@Dataset_summary$condition_levels))

    for (k in 1:length(object@Dataset_summary$condition_levels)){
      sel_mat[[k]] <- lapply(object@Stable.Networks, function(r) r$mat[[k]])
      sel_mat[[k]] <- Reduce("+", sel_mat[[k]])
      if (UnevenGroups){
        cat(paste0("Selection probabilities with subsampling ...",object@Dataset_summary$condition_levels[[k]],"...\n"))
        sel_mat[[k]] <- sel_mat[[k]]/(nCores * nreps)
      } else {
        cat(paste0("Selection probabilities without subsampling ...",object@Dataset_summary$condition_levels[[k]],"...\n"))
        sel_mat[[k]] <- sel_mat[[k]]/(2 * nCores * nreps)
      }
    }


    ###*********************************************###
    ## Estimate the partial correlation matrix
    ###*********************************************###
    n <- nrow(NormalExpression(object))
    p <- ncol(NormalExpression(object))

    Ip <- diag(rep(1,p))

    ## Model selection is done via adjusted DGlasso, where the inverse frequency weighted graphical lasso is applied.
    wAdj <- vector("list", length(object@Dataset_summary$condition_levels))
    #is necessary?
    Qmat <- vector("list", length(object@Dataset_summary$condition_levels))
    pCorMat <- vector("list", length(object@Dataset_summary$condition_levels))

    for (k in 1:length(object@Dataset_summary$condition_levels)){
      cat(paste0('Estimating model ...', object@Dataset_summary$condition_levels[k], '...\n'))
      fit <- adjDGlasso_minimal(t(object@Dataset_summary$scaled_separated_conditions[[k]]), weights=1/(1e-04 + sel_mat[[k]]))
      wAdj[[k]] <- fit$Theta.glasso
    }

    ## Get the unweighted adjacency matrix by thresholding the partial correlations
    Ahat <- NULL
    for (k in 1:length(object@Dataset_summary$condition_levels)){
      Ahat[[k]] <- abs(wAdj[[k]]) >= matrix(rep(eps, p^2), p, p)
    }

    cat(paste0("Number of edges in ",object@Dataset_summary$condition_levels[[1]],": ", sum(Ahat[[1]])/2, "\n"))
    cat(paste0("Number of edges in ",object@Dataset_summary$condition_levels[[2]],": ", sum(Ahat[[2]])/2, "\n"))

    output <- vector('list', 2)
    names(output) <-c("weighted_matrices", "threshold_matrices")
    output[["weighted_matrices"]] <- wAdj
    output[["threshold_matrices"]] <- Ahat

    object@Adjacency.Matrices <- output

    ###*********************************************###
    ##    			Output edge list              ###
    ###*********************************************###
    pairs <- combn(as.character(object@Metadata$Features), 2, simplify=FALSE)
    df <- data.frame(Metabolite.A=rep(0,length(pairs)), Metabolite.B=rep(0,length(pairs)),
                     pcor.0=rep(0,length(pairs)), pcor.1=rep(0,length(pairs)),
                     qval.0=rep(0,length(pairs)), qval.1=rep(0,length(pairs)),
                     check.names = FALSE)
    df[,1:2] <- do.call(rbind, pairs)
    df[,3] <- lowerTriangle(wAdj[[1]])
    df[,4] <- lowerTriangle(wAdj[[2]])
    df$edge <- rep(-99, length(pairs))#non-edge
    df$edge[which((abs(df$pcor.0) >= eps)*(abs(df$pcor.1) >= eps)==1)] <- "Both" #common edge
    df$edge[which((abs(df$pcor.0) >= eps)*(abs(df$pcor.1)  < eps)==1)] <- object@Dataset_summary$condition_levels[[1]]
    df$edge[which((abs(df$pcor.0) <  eps)*(abs(df$pcor.1) >= eps)==1)] <- object@Dataset_summary$condition_levels[[2]]
    df <- df[(df$edge!=-99),]
    rownames(df) <- NULL

    object@Edges <- df
    return(object)

}

runConsensusCluster <- function(object, tau0 = 0.05){
  ## Joint the two networks
  myGraph <- vector("list", length(object@Adjacency.Matrices[["threshold_matrices"]]))
  for (loop_el in 1:length(object@Adjacency.Matrices[["threshold_matrices"]])) {
    g <- graph_from_adjacency_matrix(object@Adjacency.Matrices[["weighted_matrices"]][[loop_el]], mode="undirected", weighted = TRUE)
    V(g)$name <- as.character(object@Metadata$Features)
    myGraph[[loop_el]] <- g
  }

  jointGraph <- igraph::union(myGraph[[1]], myGraph[[2]])
  jointLayout <- layout_nicely(jointGraph)
  E(jointGraph)$lty <- 1
  E(jointGraph)$color <- "black"
  E(jointGraph)$lty[is.na(E(jointGraph)$weight_2)] <- 2 #Group_1
  E(jointGraph)$lty[is.na(E(jointGraph)$weight_1)] <- 3 #Group_2
  E(jointGraph)$color[is.na(E(jointGraph)$weight_2)] <- "green" #Group_1
  E(jointGraph)$color[is.na(E(jointGraph)$weight_1)] <- "red"   #Group_2
  V(jointGraph)$color <- ifelse(object@Nodes$DEstatus=="TRUE", "purple", "white")
  V(jointGraph)$DE <- object@Nodes$DEstatus

  ###*********************************************###
  ###        Ensemble community detection 		###
  ###         with consensus clustering			###
  ###*********************************************###
  p <- numFeatures(object)
  fit <- run_consensus_cluster(jointGraph,tau=tau0,method="ensemble")
  consensus_membership <- fit$dcl
  B <- matrix(0, nrow=length(unique(consensus_membership)), p)
  rownames(B) <- paste0("Subnetwork",1:length(unique(consensus_membership)))
  for (j in 1:nrow(B)){
    B[j,which(consensus_membership==j)] <- 1
  }
  if (length(which(rowSums(B)<5))>0){
    B <- B[-which(rowSums(B)<5),]
  }
  npath <- nrow(B)

  summary_list <- list()
  for (loop_cluster in 1:nrow(B) ){
    cluster_c <- induced.subgraph(jointGraph, V(jointGraph)$name[(B[loop_cluster,]==1)])
    summary_list[[loop_cluster]] <- data.frame("number.of.nodes"=length(V(cluster_c)),
                                               "number.of.edges"=length(E(cluster_c)),
                                               "number.of.DE.nodes"=sum(as.numeric(table(V(cluster_c)$DE)[names(table(V(cluster_c)$DE))==TRUE])),
                                               "number.of.DE.edges"=sum(as.numeric(table(E(cluster_c)$color)[names(table(E(cluster_c)$color)) %in% c("red", "green")])),check.names = FALSE)
  }

  summary_stat <- data.frame("Subnetworks"= rownames(B), do.call(rbind, summary_list), check.names = FALSE)

  object@NetGSA[["summary"]]<- summary_stat
  object@NetGSA[["B"]] <- B
  object@Nodes$membership <- consensus_membership
  object@Joint.Graph <- jointGraph

  return(object)

}



runNetGSA <- function(object){

  out.netgsa <- NetGSA(object@Adjacency.Matrices[["weighted_matrices"]],
                       x = cbind(object@Dataset_summary$scaled_separated_conditions[[1]], object@Dataset_summary$scaled_separated_conditions[[2]]),
                       y = c(rep(1, ncol(object@Dataset_summary$scaled_separated_conditions[[1]])), rep(2, ncol(object@Dataset_summary$scaled_separated_conditions[[2]]))),
                       B = object@NetGSA$B, lklMethod = "REML")

  #add netGSA results to Node list
  object@Nodes$mean1 <- out.netgsa$beta[[1]]
  object@Nodes$mean2 <- out.netgsa$beta[[2]]
  object@Nodes$meanchange <- out.netgsa$beta[[2]] - out.netgsa$beta[[1]]
  object@Nodes$mc.notes <- paste(object@Dataset_summary$condition_levels[[2]], 'over', object@Dataset_summary$condition_levels[[1]])

  res <- data.frame(object@NetGSA[['summary']],
                    "NetGSA-pval"=out.netgsa$p.value,
                    "NetGSA.pFDR"=p.adjust(out.netgsa$p.value, "BH"),
                    check.names = FALSE)

  res <- res[order(res$NetGSA.pFDR),]
  rownames(res) <- 1:nrow(res)
  object@Nodes$membership[!(object@Nodes$membership %in% gsub('Subnetwork','',res$Subnetworks))] <- NA
  object@Nodes$membership <- rownames(res)[match(object@Nodes$membership, as.numeric(gsub('Subnetwork','',res$Subnetworks)))]
  object@Nodes$membership <- as.numeric(object@Nodes$membership)
  res$Subnetworks <- paste0("Subnetwork ",rownames(res))

  object@NetGSA <- object@NetGSA[object@NetGSA != 'B']
  object@NetGSA[['summary']] <- NULL
  object@NetGSA[['summary']] <- res

  return(object)
}













