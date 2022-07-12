# Differential Network Enrichment Analysis
# (DNEA) algorithm.
# Reference for algorithm and code:
# "Ma, Jing, et al. "Differential network enrichment analysis reveals novel lipid pathways in chronic kidney disease." Bioinformatics 35.18 (2019): 3441-3452."

rm(list=ls())

## Loading packages
require(gdata)
require(zoo)
require(igraph)
require(glasso)
library(glmnet)
library(corpcor)
library(dplyr)
library(foreach)
library(parallel)
library(doParallel)

source("lib2/preprocess_lib.r")
source("lib2/JSEM.R")
source("lib2/netgsa_complex.R")


main.seed <- 101
set.seed(main.seed)

## parameter setting for each of the 3 steps
## in the analysis pipeline. When performing
## a step, set only the specific parameter to
## to TRUE and ensure other steps are set to
## FALSE

#tuning parameter selection with BIC
#set to TRUE in Step 1
BIC <- FALSE

#stability selection
#set to TRUE in Step 2
SS <- FALSE
# whether stability selection
# has to be performed coupled
# with additional subsampling
# recommended for highly unbalanced
# sample group design
performSubsamling <- FALSE

# estimate partial correlation matrix (getGraph)
# and perform NetGSA (runNetGSA)
# set to TRUE in Step 3
runNetGSA <- TRUE
getGraph <- TRUE
savePlots <- TRUE


## other parameter settings

# FDR threshold for defining whether
# a compound is differential
fdr.cutoff <- 0.05

# number of cores for stability selection (Step 2)
nCores <- 25

# number of stability selection reps per core
nreps <- 20

# threshold for consensus matrix (can remain default)
tau0 <- 0.5

# number of experimental conditions in the data
ncond <- 2

eps <- 1e-06

## input data information

# sample group designation in input data
group1 <- "COVID Negative"
group2 <- "COVID Positive"

# input data filename and directory
# (specify complete path to directory)
inFolder <- "~/Documents/MS_UMICH/Karnovsky_lab/Filigree/DNEA-script/"
filename <- "COVID_data"

# output directory name
# (specify complete path to directory)
OutFolder <- "~/Documents/MS_UMICH/Karnovsky_lab/Filigree/DNEA_code_results/COVID_test_results/"


if (!getGraph){
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

  if (BIC){ ## select the tuning parameter
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

  } else if (!BIC && SS){


  }
}

  if (SS){##stability selection, which requires lastar.guo from the previous step

    listX = lapply(dat, t)

	if (performSubsamling){
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

if (getGraph){
  load(paste0(OutFolder,filename,"_stable_networks",".rda"))
  load(paste0(OutFolder,filename,"_dataset_summary",".rda"))

  ## Retrieve stable networks, which requires the stability selection results from previous step
  sel_mat <- vector("list", ncond)

  for (k in 1:ncond){
    sel_mat[[k]] <- lapply(stab_guo, function(r) r$mat[[k]])
    sel_mat[[k]] <- Reduce("+", sel_mat[[k]])
	if (performSubsamling){
		cat("Selection probabilities with subsampling ... \n")
		sel_mat[[k]] <- sel_mat[[k]]/(nCores * nreps)
	} else {
		cat("Selection probabilities without subsampling ... \n")
		sel_mat[[k]] <- sel_mat[[k]]/(2 * nCores * nreps)
	}
  }


  ###*********************************************###
  ## Estimate the partial correlation matrix
  ###*********************************************###
  n <- ncol(dataset$dat)
  p <- nrow(dataset$dat)
  x <- dataset$dat
  xx <- vector("list", ncond)
  xx[[1]] = x[, which(dataset$sample_info$Sample_Group == group1)]
  xx[[2]] = x[, which(dataset$sample_info$Sample_Group == group2)]

  Ip <- diag(rep(1,p))

  ## Model selection is done via adjusted DGlasso, where the inverse frequency weighted graphical lasso is applied.
  wAdj <- vector("list", ncond)
  Qmat <- vector("list", ncond)
  pCorMat <- vector("list", ncond)
  print('test passed')
  for (k in 1:ncond){
    cat('Estimating model ...', k, '...\n')
	fit <- adjDGlasso_minimal(t(xx[[k]]), weights=1/(1e-04 + sel_mat[[k]]))
	wAdj[[k]] <- fit$Theta.glasso
   }
  print('test passed')
  ## Get the unweighted adjacency matrix by thresholding the partial correlations
  Ahat <- NULL
  for (k in 1:ncond){
    Ahat[[k]] <- abs(wAdj[[k]]) >= matrix(rep(eps, p^2), p, p)
  }

  cat("Number of edges in Group_1: ", sum(Ahat[[1]])/2, "\n")
  cat("Number of edges in Group_2: ", sum(Ahat[[2]])/2, "\n")

  save(wAdj, Ahat, file=paste0(OutFolder,filename ,"_adjacency_matrices",".rda"))

  ###*********************************************###
  ##    			Output edge list              ###
  ###*********************************************###
  pairs <- combn(as.character(dataset$metab_info$Compound), 2, simplify=FALSE)
  df <- data.frame(Metabolite.A=rep(0,length(pairs)), Metabolite.B=rep(0,length(pairs)),
                   pcor.0=rep(0,length(pairs)), pcor.1=rep(0,length(pairs)),
                   qval.0=rep(0,length(pairs)), qval.1=rep(0,length(pairs)),
				   check.names = FALSE)
  df[,1:2] <- do.call(rbind, pairs)
  df[,3] <- lowerTriangle(wAdj[[1]])
  df[,4] <- lowerTriangle(wAdj[[2]])
  df$edge <- rep(-99, length(pairs))#non-edge
  df$edge[which((abs(df$pcor.0) >= eps)*(abs(df$pcor.1) >= eps)==1)] <- "Both" #common edge
  df$edge[which((abs(df$pcor.0) >= eps)*(abs(df$pcor.1)  < eps)==1)] <- "Group_1"
  df$edge[which((abs(df$pcor.0) <  eps)*(abs(df$pcor.1) >= eps)==1)] <- "Group_2"
  df <- df[(df$edge!=-99),]
  rownames(df) <- NULL

  write.table(df, file=paste0(OutFolder,filename,"_edgelist.txt"), row.names=FALSE, sep = "\t", quote = FALSE)

  ## Joint the two networks
  myGraph <- vector("list", length(Ahat))
  for (loop_el in 1:length(Ahat)) {
    g <- graph_from_adjacency_matrix(wAdj[[loop_el]], mode="undirected", weighted = TRUE)
    V(g)$name <- as.character(dataset$metab_info$ShortName)
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
  V(jointGraph)$color <- ifelse(dataset$metab_info$DEstatus=="TRUE", "purple", "white")
  V(jointGraph)$DE <- dataset$metab_info$DEstatus

  save(jointGraph, file=paste0(OutFolder,filename,"_joint_graph",".rda"))

  if (runNetGSA){

    ###*********************************************###
    ###        Ensemble community detection 		###
	###         with consensus clustering			###
    ###*********************************************###
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

    dataset$metab_info$membership <- consensus_membership

    ###*********************************************###
    ###                   NetGSA					###
    ###*********************************************###
    out.netgsa <- NetGSA(wAdj, x = cbind(xx[[1]], xx[[2]]), y = c(rep(1, ncol(xx[[1]])), rep(2, ncol(xx[[2]]))), B = B, lklMethod = "REML")


    # ###*********************************************###
    # ##                   Run GSA
    # ###*********************************************###
    # colnames(B) <- paste("g", 1:p, sep = "")
    # genesets <- list()
    # for (gs in 1:nrow(B)){
    #   genesets[[gs]] <- names(which(B[gs,]==1))
    # }
    # geneset.names=gsub('Cluster','set',rownames(B))
    #
    # xt = cbind(xx[[1]], xx[[2]])
    # yt = c(rep(1, ncol(xx[[1]])), rep(2, ncol(xx[[2]])))
    #
    # GSA.obj<-GSA(xt, yt, genenames=colnames(B), genesets=genesets, resp.type="Two class unpaired", nperms=3000,restand = F,minsize = 5,s0=1)
    # GSAouts <- GSA.listsets(GSA.obj, geneset.names = geneset.names, FDRcut = 1)
    # GSAouts = rbind(GSAouts$negative[, c(1:5)], GSAouts$positive[, c(1:5)])
    #
    # nsets = length(genesets)
    # indx = as.numeric(c(1:nsets)[!(geneset.names %in% GSAouts[, 2])])
    # toadd = cbind(indx, geneset.names[indx], rep(0, length(indx)), rep(1, length(indx)), rep(1, length(indx)))
    # if (length(indx) > 0) {
    #   GSAouts = rbind(GSAouts, toadd)
    # }
    # ord = order(as.numeric(GSAouts[, 1]))
    # GSAoutsF = GSAouts[ord, ]
    #
    # #with restandardization
    # GSA.obj<-GSA(xt, yt, genenames=colnames(B), genesets=genesets, resp.type="Two class unpaired", nperms=3000,restand = T,minsize = 5,s0=1)
    # GSAouts <- GSA.listsets(GSA.obj, geneset.names = geneset.names, FDRcut = 1)
    # GSAouts = rbind(GSAouts$negative[, c(1:5)], GSAouts$positive[, c(1:5)])
    #
    # nsets = length(genesets)
    # indx = as.numeric(c(1:nsets)[!(geneset.names %in% GSAouts[, 2])])
    # toadd = cbind(indx, geneset.names[indx], rep(0, length(indx)), rep(1, length(indx)), rep(1, length(indx)))
    # if (length(indx) > 0) {
    #   GSAouts = rbind(GSAouts, toadd)
    # }
    # ord = order(as.numeric(GSAouts[, 1]))
    # GSAoutsT = GSAouts[ord, ]
    #


    ## Output node information
	dataset$metab_info$mean1 <- out.netgsa$beta[[1]]
	dataset$metab_info$mean2 <- out.netgsa$beta[[2]]
    dataset$metab_info$meanchange <- out.netgsa$beta[[2]] - out.netgsa$beta[[1]]
	dataset$metab_info$mc.notes <- "Group_2 over Group_1"

    ###*********************************************###
    ##        Output enrichment results
    ###*********************************************###
    res <- data.frame(summary_stat,
                      "NetGSA-pval"=out.netgsa$p.value, "NetGSA.pFDR"=p.adjust(out.netgsa$p.value, "BH"), check.names = FALSE)
                      # "GSA-F-pval"=as.numeric(GSAoutsF[, 4]), "GSA-F-pFDR"=as.numeric(GSAoutsF[, 5]),
                      # "GSA-T-pval"=as.numeric(GSAoutsF[, 4]), "GSA-T-pFDR"=as.numeric(GSAoutsT[, 5]))

	res <- res[order(res$NetGSA.pFDR),]
    rownames(res) <- 1:nrow(res)
    dataset$metab_info$membership <- consensus_membership
    dataset$metab_info$membership[!(dataset$metab_info$membership %in% gsub('Subnetwork','',res$Subnetworks))] <- NA
    dataset$metab_info$membership <- rownames(res)[match(dataset$metab_info$membership, as.numeric(gsub('Subnetwork','',res$Subnetworks)))]
    dataset$metab_info$membership <- as.numeric(dataset$metab_info$membership)
    res$Subnetworks <- paste0("Subnetwork ",rownames(res))

    write.csv(res, file=paste0(OutFolder,filename ,"_netgsa.csv"), row.names=FALSE)

    write.csv(dataset$metab_info,
	file=paste0(OutFolder,filename ,"_nodelist.csv"),
              row.names=FALSE)

	save.image(file=paste0(OutFolder,filename,"_netgsa_results.rda"))

    ## Output the clusters
    if (savePlots){
      ## Visualize the clusters and the mean changes
      V(jointGraph)$size <- abs(dataset$metab_info$meanchange)*30
      pdf(paste0(OutFolder,filename,"_consensus_clusters_fdr",fdr.cutoff*100,'_',".pdf"), height = 10, width = 10)
      for (loop_cluster in 1:nrow(B) ){
        cluster_c <- induced.subgraph(jointGraph, V(jointGraph)$name[which(dataset$metab_info$membership==loop_cluster)])
        plot(cluster_c, vertex.label = V(cluster_c)$name, vertex.label.cex = 1,
             layout = layout.fruchterman.reingold(cluster_c),
             main = paste0(res$cluster.name[loop_cluster]," ( qvalue - ",round(res$NetGSA.pFDR[loop_cluster],2),")"))
      }
      dev.off()
    }
  }
}
