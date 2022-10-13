#' @include JSEM.R
#'
NULL
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Optimize the Lambda parameter for glasso
#'
#' This function will calculate the Bayesian information criterion (BIC) and liklihood for a range of lambda values
#'  determined by the number of features within the dataset. The lambda value with the lowest (BIC) score is
#'  chosen for analysis. It takes a DNEAobject as input and has an option that allows the function to be run
#'  in parallel. The main.seed parameter is for reproducibility and can be left as default.
#'
#' @param object A DNEA object
#' @param nCores The number of cores available for parallel processing. If more cores than lambda values
#'        tested is specified, will default to one worker per lambda value.If set to 1, the analysis is not
#'        run in parallel.
#'@param lambda_values An optional list of lambda values to fit a model and calculate the BIC score.
#'       If not provided, a set of lambda values are chosen based on the size of the dataset.
#'@param eps_cutoff A significance cut-off for thresholding network interactions.
#'       The default value is 1e-06.
#'@param eta_value default parameter ??. Default is 0.1
#'
#' @return A DNEAobject containing the BIC and liklihood scores for every lambda value tested, as well as
#'         the optimized lambda value
#' @export
#' @import zoo
#' @import glasso
#' @import glmnet
#' @import corpcor
#' @import dplyr
#' @import parallel
#' @import pbapply
BICtune <- function(object,
                    nCores = 1,
                    lambda_values,
                    eps_cutoff = 1e-06,
                    eta_value = 0.1){

  ############################################
  #**Prepare data and initialize parameters**#
  ############################################

  #split data
  dat <- split_by_condition(dat = scaledExpressionData(object),
                            condition_levels = object@dataset_summary[['condition_levels']],
                            condition_by_sample = condition(object))

  ## Joint estimation
  n4cov <- max(sapply(dat, ncol))

  #Pre-define a range of lambda to select the tuning parameters if none are provided
  if(missing(lambda_values)){
    lambda_values <- seq(0.01, 0.3, 0.02)*sqrt(log(numFeatures(object))/n4cov)
  }else{
    lambda_values <- unlist(lambda_values)
  }

  #This needs to be more informative based on the data.
  trainX <- t(do.call(cbind, dat))
  trainY <- c(rep(1, ncol(dat[[1]])),
              rep(2, ncol(dat[[2]])))

  #############################################
  #**Initialize workers and optimize lambda **#
  #############################################

  message("BIC using Guo et al ...", appendLF = TRUE)

  #If more cores available than lambda values tested, will set one worker
  #per lambda value to improve efficiency
  if(nCores > 1){

    message("Lambda parameter will be optimized in parallel.\n", appendLF = TRUE)

    if(nCores > length(lambda_values)){

      nCores <- length(lambda_values)
      message("More cores available than lambda values being tested.")
      message("Resource utilization is being adjusted for efficiency", appendLF = TRUE)
    }

    #initialize parallel process
    cl <- parallel::makeCluster(nCores)
    on.exit(stopCluster(cl))

    #set progress bar
    pbapply::pboptions(type = "timer", char = c('='), style = 5)

    #pass necessary objects to workers
    parallel::clusterExport(cl = cl,
                            varlist = c("lambda_values",
                                        "trainX",
                                        "trainY",
                                        "eta_value"),
                            envir = environment())
    parallel::clusterEvalQ(cl = cl, c(library("MASS"), library("glasso")))
    parallel::clusterExport(cl = cl, varlist = c("CGM_AHP_tune","CGM_AHP_train", "matTr"))

    # optimize lambda
    bic_guo <- pblapply(cl = cl,
                        FUN = 'CGM_AHP_tune',
                        X = lambda_values,
                        trainX = trainX,
                        testX = trainX,
                        model = trainY,
                        BIC = TRUE,
                        eps = eps_cutoff,
                        eta = eta_value)

  } else{

    message('Lambda parameter will be optimized sequentially. \n', appendLF = TRUE)

    # optimize lambda
    bic_guo <- pblapply(FUN = 'CGM_AHP_tune',
                        X = lambda_values,
                        trainX = trainX,
                        testX = trainX,
                        model = trainY,
                        BIC = TRUE,
                        eps = eps_cutoff,
                        eta = eta_value)

  }

  ##make sure all bic values are finite and remove those that are not
  tmp = unlist(sapply(bic_guo, function(a) a$BIC))
  if (max(is.infinite(tmp))==1){
    bic_guo <- bic_guo[is.finite(tmp)]
    lambda_values <- lambda_values[is.finite(tmp)]
    lastar_guo <- lambda_values[which.min(tmp)]
  } else {
    lastar_guo <- lambda_values[which.min(sapply(bic_guo, function(a) a$BIC))]
  }

  #######################
  #**update DNEAobject**#
  #######################

  object@hyperparameter[['BIC_scores']] <- bic_guo
  object@hyperparameter[['optimized_lambda']] <- lastar_guo
  object@hyperparameter[['tested_lambda_values']] <- lambda_values

  message(paste0('The optimal Lambda hyper-parameter has been set to: ', lastar_guo, '!'))
  return(object)
}

#' Performs stability selection to determine probability of feature interactions
#'
#' This function randomly samples each condition a number of times specified by nreps. It takes a DNEAobject as input.
#' The unevenGroups parameter should be set to TRUE if the sample numbers across condition are uneven - this will
#' ensure that samples are randomly sampled in a way that ensures equal representation in stability selection. This
#' function can also be run in parallel by specifying the available cores through nCores. The main.seed parameter
#' is for reproducibility and can be left as default.
#'
#' @param object A DNEA object
#' @param runParallel A boolean indicating if stability selection should be run in parallel
#' @param subSample A boolean that specifies whether the number of samples are unevenly split
#'         by condition and, therefore, should be adjusted for when randomly sampling.
#' @param nreps The total number of reps to perform stability selection. As the sample number
#'        gets smaller, more reps should be run; For datasets under 200 samples, we suggest running
#'        500 reps.
#' @param nCores The number of cores available for parallel processing. It is generally optimal
#'        to begin with ~10 reps per core. However, you may need to run one rep per core for optimal
#'        performance with some datasets.
#' @param optimized_lambda The optimal lambda value to be used in the model. This parameter is only
#'        necessary if BICtune() is not performed
#' @param main.seed Sets the seed for random number generation. This ensures reproducibility.
#'
#' @return DNEAobject containing the stable networks for analysis.
#'
#' @export
#' @import zoo
#' @import stats
#' @import glasso
#' @import glmnet
#' @import corpcor
#' @import dplyr
#' @import parallel
#' @import pbapply
stabilitySelection <- function(object,
                               runParallel = TRUE,
                               subSample = FALSE,
                               nreps = 50,
                               nCores = 4,
                               optimal_lambda,
                               main.seed = 101){

  ############################################
  #**Prepare data and initialize parameters**#
  ############################################

  #stability selection requires lambda hyper-parameter. Will use optimal_lambda if
  #supplied, otherwise looks for @hyperparameter[["optimized_lambda"]] in DNEAobject

  if(is.null(object@hyperparameter[['optimized_lambda']]) == FALSE){
    if(missing(optimal_lambda) == FALSE){

      optimized_lambda <- optimal_lambda
      warning('optimal_lambda argument was provided even though @hyperparameter[["optimized_lambda"]]
              already exists - optimal_lambda will be used in analysis')

    } else{

      optimized_lambda <- object@hyperparameter[["optimized_lambda"]]
    }
  } else{
    if(missing(optimal_lambda) == FALSE){

      optimized_lambda <- optimal_lambda
      object@hyperparameter[["optimized_lambda"]] <- optimal_lambda
      message('@hyperparameter[["optimized_lambda"]] was previously empty and now set to optimal_lambda')

    } else{

      stop('No lambda value was supplied for the model. Please run BICtune() or provide a lambda
         value using the optimal_lambda parameter.')
    }
  }

  #split data by condition
  data_split_by_condition = lapply(split_by_condition(dat = scaledExpressionData(object),
                                    condition_levels = object@dataset_summary[['condition_levels']],
                                    condition_by_sample = condition(object)), function(d) t(d))


  # pbapply requires a vector of length nreps
  nreps_input = 1:nreps

  #initialize static variables to pass to workers
  stabsel_init_param <- stabsel_init(listX = data_split_by_condition, nreps = nreps)

  #########################################################
  #**Initialize workers and perform stability selection **#
  #########################################################

  #set progress bar
  pbapply::pboptions(type = "timer", char = c('='), style = 5)

  #print lambda used
  message(paste0('Using Lambda hyper-parameter: ', optimized_lambda,'!'))

  if(runParallel){

    message('stabilitySelection() will be run in parallel ...', appendLF = TRUE)

    #create independent processes to run reps in parallel
    cl <- makeCluster(nCores)
    on.exit(stopCluster(cl))

    #pass necessary objects to independent workers
    parallel::clusterExport(cl = cl, varlist = c("CGM_AHP_tune","CGM_AHP_train","CGM_AHP_stabsel_subsample","CGM_AHP_stabsel", "matTr"))
    parallel::clusterExport(cl = cl, varlist = c("data_split_by_condition", "optimized_lambda","nreps","main.seed","stabsel_init_param"), envir = environment())
    parallel::clusterEvalQ(cl = cl, c(library("MASS"), library("glasso"), set.seed(main.seed)))


    if (subSample){
      message("Stability selection WITH additional subsampling using Guo et al ...\n", appendLF = TRUE)

      stab_sel <- pbapply::pblapply(cl = cl,
                                    FUN = "CGM_AHP_stabsel_subsample",
                                    X = nreps_input,
                                    init_param = stabsel_init_param,
                                    listX = data_split_by_condition,
                                    lastar = optimized_lambda)
    } else {
      message("Stability selection WITHOUT additional subsampling using Guo et al ...\n", appendLF = TRUE)

      stab_sel <- pbapply::pblapply(cl = cl,
                                    FUN ="CGM_AHP_stabsel",
                                    X = nreps_input,
                                    init_param = stabsel_init_param,
                                    listX = data_split_by_condition,
                                    lastar = optimized_lambda)
    }
  }else{

    message('stabilitySelection() will be run sequentially ...', appendLF = TRUE)

    if (subSample){
      message("Stability selection WITH additional subsampling using Guo et al ...\n", appendLF = TRUE)

      stab_sel <- pbapply::pblapply(FUN = "CGM_AHP_stabsel_subsample",
                                    X = nreps_input,
                                    init_param = stabsel_init_param,
                                    listX = data_split_by_condition,
                                    lastar = optimized_lambda)
    } else {
      message("Stability selection WITHOUT additional subsampling using Guo et al ...\n", appendLF = TRUE)

      stab_sel <- pbapply::pblapply(FUN ="CGM_AHP_stabsel",
                                    X = nreps_input,
                                    init_param = stabsel_init_param,
                                    listX = data_split_by_condition,
                                    lastar = optimized_lambda)
    }
  }

  #####################################
  #**Concatenate results for output **#
  #####################################

  #initiate list for stability selection raw results
  selection_results <- vector("list", length(object@dataset_summary$condition_levels))
  names(selection_results) <- object@dataset_summary$condition_levels

  #initiate list for stability selection results converted to probabilities
  selection_probabilities <- vector("list", length(object@dataset_summary$condition_levels))
  names(selection_probabilities) <- object@dataset_summary$condition_levels

  #initiate output list
  output_stabsel <- vector("list", length = 2)
  names(output_stabsel) <- c("selectionResults","selectionProbabilites")

  #reduce results to one matrix and calculate selection probabilities
  for (k in 1:length(object@dataset_summary$condition_levels)){
    selection_results[[k]] <- lapply(stab_sel, function(r) r$mat[[k]])
    selection_results[[k]] <- Reduce("+", selection_results[[k]])
    if (subSample){
      message(paste0("Calculating selection probabilities WITH subsampling for...",object@dataset_summary$condition_levels[[k]],"..."), appendLF = TRUE)
      selection_probabilities[[k]] <- selection_results[[k]]/(nreps)
    } else {
      message(paste0("Calculating selection probabilities WITHOUT subsampling for...",object@dataset_summary$condition_levels[[k]],"..."), appendLF = TRUE)
      selection_probabilities[[k]] <- selection_results[[k]]/(2 * nreps)
    }
  }

  output_stabsel[["selectionResults"]] <- selection_results
  output_stabsel[["selectionProbabilites"]] <- selection_probabilities

  object@stable_networks <- output_stabsel
  return(object)
}

#' Creates the network model using glasso
#'
#' This function takes in a DNEAobjct and fits a glasso model using the optimized lambda value
#' determined via BICtune() or otherwise specified. If selection probabilites for each node were
#' calculated using stabilitySelection(), those are also utilized. An adjacency matrix and a
#' thresholded adjacency matrix, determiend by the eps value, is output.
#'
#' @param object A DNEAobject
#' @param eps A significance cut-off for thresholding the adjacency matrix. This can be left as default
#'
#' @return A DNEA object containing an adjacency matrix for the data network,
#'          determined by the glasso model
#'
#' @import zoo
#' @import gdata
#' @import glasso
#' @import glmnet
#' @import corpcor
#' @import dplyr
#' @export
getNeworks <- function(object, eps = 1e-06){

  ############################################
  #**Prepare data and initialize parameters**#
  ############################################

    num_samples <- nrow(scaledExpressionData(object))
    num_features <- ncol(scaledExpressionData(object))

    Ip <- diag(rep(1,num_features))

    ## Model selection is done via adjusted DGlasso, where the inverse frequency weighted graphical lasso is applied.
    weights_adjusted <- vector("list", length(object@dataset_summary$condition_levels))

    #separate the data by condition
    data_split_by_condition <- split_by_condition(dat = scaledExpressionData(object),
                                                    condition_levels = object@dataset_summary$condition_levels,
                                                    condition_by_sample = object@metadata$condition_values)

    #############################################
    #**Estimate the partial correlation matrix**#
    #############################################

    for (k in 1:length(object@dataset_summary$condition_levels)){
      message(paste0('Estimating model for ', object@dataset_summary$condition_levels[k], ' ...'), appendLF = TRUE)
      fit <- adjDGlasso_minimal(t(data_split_by_condition[[k]]), weights=1/(1e-04 + object@stable_networks$selectionProbabilites[[k]]))
      weights_adjusted[[k]] <- fit$Theta.glasso
    }

    ## Get the unweighted adjacency matrix by thresholding the partial correlations
    weights_threshold <- NULL
    for (k in 1:length(object@dataset_summary$condition_levels)){
      weights_threshold[[k]] <- abs(weights_adjusted[[k]]) >= matrix(rep(eps, num_features^2), num_features, num_features)
    }

    #####################################
    #**Concatenate results for output **#
    #####################################

    #print message for total edges
    message(paste0("Number of edges in ",object@dataset_summary$condition_levels[[1]],": ", sum(weights_threshold[[1]])/2), appendLF = TRUE)
    message(paste0("Number of edges in ",object@dataset_summary$condition_levels[[2]],": ", sum(weights_threshold[[2]])/2), appendLF = TRUE)

    output <- vector('list', 2)
    names(output) <-c("weighted_matrix", "threshold_matrix")
    output[["weighted_matrix"]] <- weights_adjusted
    output[["threshold_matrix"]] <- weights_threshold

    object@adjacency_matrix <- output

    #######################
    #**Create Edge List **#
    #######################

    #initiate output dataframe
    pairs <- combn(as.character(object@metadata$features), 2, simplify=FALSE)
    df <- data.frame(Metabolite.A=rep(0,length(pairs)), Metabolite.B=rep(0,length(pairs)),
                     pcor.0=rep(0,length(pairs)), pcor.1=rep(0,length(pairs)),
                     check.names = FALSE)

    #concatenate results into dataframe
    df[,1:2] <- do.call(rbind, pairs)
    df[,3] <- lowerTriangle(weights_adjusted[[1]])
    df[,4] <- lowerTriangle(weights_adjusted[[2]])
    df$edge <- rep(-99, length(pairs))#non-edge
    df$edge[which((abs(df$pcor.0) >= eps)*(abs(df$pcor.1) >= eps)==1)] <- "Both" #common edge
    df$edge[which((abs(df$pcor.0) >= eps)*(abs(df$pcor.1)  < eps)==1)] <- object@dataset_summary$condition_levels[[1]]
    df$edge[which((abs(df$pcor.0) <  eps)*(abs(df$pcor.1) >= eps)==1)] <- object@dataset_summary$condition_levels[[2]]
    df <- df[(df$edge!=-99),]
    rownames(df) <- NULL

    object@edge_list <- df
    return(object)

}
#' runConsensusCluster performs consensus clustering using an adjacency matrix for the network
#'
#' This function will take as input an adjacency matrix from the determined networks and perform
#' consensus clustering using the following methods from the igraph package: cluster_edge_betweenness,
#' cluster_fast_greedy, cluster_infomap, cluster_label_prop, cluster_leading_eigen, cluster_louvain,
#' cluster_walktrap. The output results in sub-network classification for the nodes within the network.
#'
#' @param object A DNEA object
#' @param tau The consensus probabilty threshold for agreement among clustering algorithms
#' @param method The consensus clustering method to be used. The options are as follows
#'        - "ensemble": indicates that all seven of the available clustering methods
#'        (cluster_edge_betweenness, cluster_fast_greedy, cluster_infomap, cluster_label_prop,
#'        cluster_leading_eigen, cluster_louvain, cluster_walktrap) should be used.
#'        -"lpm" utilizes the cluster_label_prop, cluster_infomap, and cluster_walktrap methods.
#'        -"walktrap" utilizes only cluster_walktrap.
#'        -"infomap" utilizes only infomap
#' @param num_iterations The number of clustering iterations to perform - this parameter not relevant
#'        for the "ensemble" method. Default is 10 iterations.
#'
#' @return A DNEAobject containing sub-network determinations for the nodes within the input network.
#'        A summary of the consensus clustering results can be viewed using getClusterResults().
#'        Sub-network classification for each node can be found in the node_list slot of the returned
#'        DNEAobject.
#'
#' @import igraph
#' @export
runConsensusCluster <- function(object, tau = 0.5, num_iterations = 10, method = "ensemble"){

  #####################################
  #**Join the two condition networks**#
  #####################################

  #create list to hold graph from adjacency matrix
  adjacency_matrix_graphs <- vector("list", length(object@adjacency_matrix[["threshold_matrix"]]))

  for (loop_el in 1:length(object@adjacency_matrix[["threshold_matrix"]])) {
    g <- graph_from_adjacency_matrix(object@adjacency_matrix[["weighted_matrix"]][[loop_el]], mode="undirected", weighted = TRUE)
    V(g)$name <- as.character(object@metadata$features)
    adjacency_matrix_graphs[[loop_el]] <- g
  }

  #join adjacency matrix graphs
  jointGraph <- igraph::union(adjacency_matrix_graphs[[1]], adjacency_matrix_graphs[[2]])

  #modify the jointGraph
  jointLayout <- layout_nicely(jointGraph)
  E(jointGraph)$lty <- 1
  E(jointGraph)$color <- "black"
  E(jointGraph)$lty[is.na(E(jointGraph)$weight_2)] <- 2 #Group_1
  E(jointGraph)$lty[is.na(E(jointGraph)$weight_1)] <- 3 #Group_2
  E(jointGraph)$color[is.na(E(jointGraph)$weight_2)] <- "green" #Group_1
  E(jointGraph)$color[is.na(E(jointGraph)$weight_1)] <- "red"   #Group_2
  V(jointGraph)$color <- ifelse(object@node_list$DEstatus=="TRUE", "purple", "white")
  V(jointGraph)$DE <- object@node_list$DEstatus

  ###########################################################
  #**ensemble community detection with consusus clustering**#
  ###########################################################

  #run consensus cluster algorithm
  # fit <- run_consensus_cluster(jointGraph,tau=tau0,method="ensemble", runParallel = runParallel, nCores = nCores)
  fit <- run_consensus_cluster(jointGraph, tau=tau, method = method, num_iterations = num_iterations )
  consensus_membership <- fit$dcl

  #gather results
  subnetwork_results <- matrix(0, nrow=length(unique(consensus_membership)), numFeatures(object))
  rownames(subnetwork_results) <- paste0("Subnetwork",1:length(unique(consensus_membership)))
  for (j in 1:nrow(subnetwork_results)){
    subnetwork_results[j,which(consensus_membership==j)] <- 1
  }
  if (length(which(rowSums(subnetwork_results)<5))>0){
    subnetwork_results <- subnetwork_results[-which(rowSums(subnetwork_results)<5),]
  }
  npath <- nrow(subnetwork_results)

  #####################################
  #**Concatenate results for output **#
  #####################################

  summary_list <- list()
  for (loop_cluster in 1:nrow(subnetwork_results) ){
    cluster_c <- induced.subgraph(jointGraph, V(jointGraph)$name[(subnetwork_results[loop_cluster,]==1)])
    summary_list[[loop_cluster]] <- data.frame("number_of_nodes"=length(V(cluster_c)),
                                               "number_of_edges"=length(E(cluster_c)),
                                               "number_of_DE.nodes"=sum(as.numeric(table(V(cluster_c)$DE)[names(table(V(cluster_c)$DE))==TRUE])),
                                               "number_of_DE.edges"=sum(as.numeric(table(E(cluster_c)$color)[names(table(E(cluster_c)$color)) %in% c("red", "green")])),check.names = FALSE)
  }

  summary_stat <- data.frame("Subnetworks"= rownames(subnetwork_results), do.call(rbind, summary_list), check.names = FALSE)

  object@netGSA_results[["summary"]]<- summary_stat
  object@netGSA_results[["subnetwork_results"]] <- subnetwork_results
  object@node_list$membership <- consensus_membership
  object@joint_graph <- jointGraph

  return(object)

}
#' Runs the NetGSA algorithm described in Hellstern et al.
#'
#' runNetGSA takes as input a DNEA object and utilizes the adjacency matrix calculated in getNetworks()
#' to perform NetGSA on the two networks to determine differences in subnetwork (pathway) expression.
#'
#' @param object A DNEAobject
#'
#' @returns A DNEAobject containing containing results from NetGSA. Pathway expression differences for
#'          each node can be found in the node_list. A summary of the NetGSA results can be viewed
#'          using getNetGSAresults().
#'
#'
#' @export
runNetGSA <- function(object){

  #################################
  #**Prepare data and run netGSA**#
  #################################
  #separate the data by condition
  separated_conditions_data <- split_by_condition(dat = scaledExpressionData(object),
                                                  condition_levels = object@dataset_summary$condition_levels,
                                                  condition_by_sample = object@metadata$condition_values)

  out.netgsa <- NetGSA(object@adjacency_matrix[["weighted_matrix"]],
                       x = cbind(separated_conditions_data[[1]], separated_conditions_data[[2]]),
                       y = c(rep(1, ncol(separated_conditions_data[[1]])), rep(2, ncol(separated_conditions_data[[2]]))),
                       B = object@netGSA_results[["subnetwork_results"]], lklMethod = "REML")

  #####################################
  #**Concatenate results for output **#
  #####################################

  #add netGSA results to Node list
  object@node_list$mean1 <- out.netgsa$beta[[1]]
  object@node_list$mean2 <- out.netgsa$beta[[2]]
  object@node_list$meanchange <- out.netgsa$beta[[2]] - out.netgsa$beta[[1]]
  object@node_list$mc.notes <- paste(object@dataset_summary$condition_levels[[2]], 'over', object@dataset_summary$condition_levels[[1]])

  #concatenate netGSA summary output
  res <- data.frame(object@netGSA_results[['summary']],
                    "NetGSA_pval"= out.netgsa$p.value,
                    "NetGSA_pFDR"= p.adjust(out.netgsa$p.value, "BH"),
                    check.names = FALSE)
  res <- res[order(res$NetGSA_pFDR),]
  rownames(res) <- 1:nrow(res)

  #modify membership column of Node list
  object@node_list$membership[!(object@node_list$membership %in% gsub('Subnetwork','',res$Subnetworks))] <- NA
  object@node_list$membership <- rownames(res)[match(object@node_list$membership, as.numeric(gsub('Subnetwork','',res$Subnetworks)))]
  object@node_list$membership <- as.numeric(object@node_list$membership)
  res$Subnetworks <- paste0("Subnetwork ",rownames(res))

  #update DNEAobject
  object@netGSA_results <- object@netGSA_results[object@netGSA_results != 'B']
  object@netGSA_results[['summary']] <- NULL
  object@netGSA_results[['summary']] <- res

  return(object)
}













