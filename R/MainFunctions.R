#' @include JSEM.R
#' @include utilities.R
#' @include preprocess_lib.R
#' @include all-methods.R
#' @include all-generics.R
#' @include all-classes.R
#'
NULL

#' BICtune Optimizes the Lambda parameter for glasso
#'
#' This function will calculate the Bayesian information criterion (BIC) and liklihood for a range of lambda values
#' determined by the number of features within the dataset. The lambda value with the lowest (BIC) score is
#' chosen for analysis. It takes a DNEAobject as input and has an option that allows the function to be run
#' in parallel.
#'
#' @param object A DNEA object
#' @param lambda_values An optional list of lambda values to fit a model and calculate the BIC score.
#'        If not provided, a set of lambda values are chosen based on the size of the dataset.
#' @param eps_threshold A significance cut-off for thresholding network edges
#'        The default value is 1e-06.
#' @param eta_value default parameter ??. Default is 0.1
#' @param BPPARAM A BiocParallel object
#'
#' @return A DNEAobject containing the BIC and liklihood scores for every lambda value tested, as well as
#'         the optimized lambda value
#'
#' @include JSEM.R
#' @include utilities.R
#' @import zoo
#' @import glasso
#' @import glmnet
#' @import corpcor
#' @import BiocParallel
#' @export
BICtune <- function(object,
                    lambda_values,
                    eps_threshold = 1e-06,
                    eta_value = 0.1,
                    BPPARAM = bpparam()){

  ############################################
  #**Prepare data and initialize parameters**#
  ############################################

  ##split data by condition
  dat <- split_by_condition(dat = expressionData(object, type = "normalized"),
                            condition_levels = networkGroups(object),
                            condition_by_sample = networkGroupIDs(object))

  ##create input for model training
  n4cov <- max(sapply(dat, ncol))
  trainX <- t(do.call(cbind, dat))
  trainY <- c(rep(1, ncol(dat[[1]])),
              rep(2, ncol(dat[[2]])))


  ##Pre-define a range of lambda values to evaluate during optimization if none are provided
  if(missing(lambda_values)){
    lambda_values <- seq(0.01, 0.3, 0.02)*sqrt(log(numFeatures(object))/n4cov)
  }else{
    lambda_values <- unlist(lambda_values)
  }
  #############################################
  #**Initialize workers and optimize lambda **#
  #############################################

  message("Optimizing the lambda hyperparameter using Bayesian-Information
          Criterion outlined in Guo et al. (2012)", appendLF = TRUE)

  BIC_guo <- BiocParallel::bplapply(X = lambda_values,
                                    FUN = 'CGM_AHP_tune',
                                    trainX = trainX,
                                    testX = trainX,
                                    model = trainY,
                                    BIC = TRUE,
                                    eps = eps_threshold,
                                    eta = eta_value,
                                    BPPARAM = BPPARAM,
                                    BPOPTIONS = bpoptions(progressbar = TRUE, tasks = 10))

  # #set progress bar
  # pbapply::pboptions(type = "timer", char = c('='), style = 5)
  #
  # #If more cores available than lambda values tested, will set one worker
  # #per lambda value to improve efficiency
  # if(runParallel){
  #
  #   message("Lambda parameter will be optimized in parallel.\n", appendLF = TRUE)
  #
  #   if(nCores > length(lambda_values)){
  #
  #     nCores <- length(lambda_values)
  #     message("More cores available than lambda values being tested.", appendLF = TRUE)
  #     message("Resource utilization is being adjusted for efficiency.", appendLF = TRUE)
  #     message(paste0(nCores, 'independent processes will be used in parallelization.'), appendLF = TRUE)
  #   }
  #
  #   #initialize parallel process
  #   cl <- parallel::makeCluster(nCores, type = 'PSOCK')
  #   on.exit(stopCluster(cl))
  #
  #   #pass necessary objects to workers
  #   parallel::clusterExport(cl = cl, varlist = c("CGM_AHP_tune","CGM_AHP_train", "matTr"))
  #   parallel::clusterEvalQ(cl = cl, c(library("MASS"), library("glasso")))
  #
  #   # optimize lambda
  #   BIC_guo <- pblapply(cl = cl,
  #                       FUN = 'CGM_AHP_tune',
  #                       X = lambda_values,
  #                       trainX = trainX,
  #                       testX = trainX,
  #                       model = trainY,
  #                       BIC = TRUE,
  #                       eps = eps_threshold,
  #                       eta = eta_value)
  #
  # } else{
  #
  #   message('Lambda parameter will be optimized sequentially. \n', appendLF = TRUE)
  #
  #   # optimize lambda
  #   BIC_guo <- pblapply(FUN = 'CGM_AHP_tune',
  #                       X = lambda_values,
  #                       trainX = trainX,
  #                       testX = trainX,
  #                       model = trainY,
  #                       BIC = TRUE,
  #                       eps = eps_threshold,
  #                       eta = eta_value)
  #
  # }

  ##add empty line after progress bar
  message("", appendLF = TRUE)

  ##collect BIC scores
  BIC_scores = unlist(sapply(BIC_guo, function(a) a$BIC))

  ##make sure all bic values are finite and remove those that are not
  if (max(is.infinite(BIC_scores)) == 1){

    BIC_guo <- BIC_guo[is.finite(BIC_scores)]
    lambda_values <- lambda_values[is.finite(BIC_scores)]
    lastar_guo <- lambda_values[match(min(BIC_scores), BIC_scores)]
  } else {

    lastar_guo <- lambda_values[match(min(BIC_scores), BIC_scores)]
  }

  #######################
  #**update DNEAobject**#
  #######################

  BICscores(object) <- BIC_guo
  optimizedLambda(object) <- lastar_guo
  lambdas2Test(object) <- lambda_values

  message(paste0('The optimal Lambda hyper-parameter has been set to: ', lastar_guo, '!'))
  return(object)
}

#' Performs stability selection to determine probability of feature interactions
#'
#' This function randomly samples each condition a number of times specified by nreps. It takes a DNEAobject as input.
#' The unevenGroups parameter should be set to TRUE if the sample numbers across condition are uneven - this will
#' ensure that samples are randomly sampled in a way that ensures equal representation in stability selection. This
#' function can also be run in parallel by specifying the available cores through nCores.
#'
#' @param object A DNEA object
#' @param subSample A boolean that specifies whether the number of samples are unevenly split
#'         by condition and, therefore, should be adjusted for when randomly sampling.
#' @param nreps The total number of reps to perform stability selection. As the sample number
#'        gets smaller, more reps should be run; For datasets under 200 samples, we suggest running
#'        500 reps.
#' @param optimal_lambda The optimal lambda value to be used in the model. This parameter is only
#'        necessary if BICtune() is not performed
#'
#' @param BPPARAM a BiocParallel object
#'
#' @return DNEAobject containing the stable networks for analysis.
#'
#' @import zoo
#' @import glasso
#' @import glmnet
#' @import corpcor
#' @import BiocParallel
#' @include JSEM.R
#' @include utilities.R
#' @export
stabilitySelection <- function(object,
                               subSample = FALSE,
                               nreps = 50,
                               optimal_lambda,
                               BPPARAM = bpparam()){

  ############################################
  #**Prepare data and initialize parameters**#
  ############################################

  #stabilitySelection requires lambda hyper-parameter. Will use optimal_lambda if
  #supplied, otherwise looks for @hyperparameter[["optimized_lambda"]] in DNEAobject
  #
  # choosing lambda follows the following algorithm:
  # 1. if @hyperparamater[['optimized_lambda']] and optimal_lambda provided, optimal_lambda is used
  #    for analysis
  # 2. if @hyperparamater[['optimized_lambda']] and optimal_lambda missing, use
  #    @hyperparamater[['optimized_lambda']]
  # 3. if @hyperparamater[['optimized_lambda']] and optimal_lambda provided, optimal_lambda is used
  #    for analysis
  # 4. if both @hyperparamater[['optimized_lambda']] and optimal_lambda are missing, throw error.
  if(!is.null(optimizedLambda(object))){
    if(!missing(optimal_lambda)){

      optimized_lambda <- optimal_lambda
      warning('optimal_lambda argument was provided even though @hyperparameter[["optimized_lambda"]]
              already exists - optimal_lambda will be used in analysis')

    } else{

      optimized_lambda <- optimizedLambda(object)
    }
  } else{
    if(!missing(optimal_lambda)){

      optimized_lambda <- optimal_lambda
      optimizedLambda(object) <- optimal_lambda
      message('@hyperparameter[["optimized_lambda"]] was previously empty and now set to the optimal_lambda provided!')

    } else{

      stop('No lambda value was supplied for the model. Please run BICtune() or provide a lambda
         value using the optimal_lambda parameter.')
    }
  }

  #split data by condition
  data_split_by_condition = lapply(split_by_condition(dat = expressionData(object, type = "normalized"),
                                    condition_levels = networkGroups(object),
                                    condition_by_sample = networkGroupIDs(object)), function(d) t(d))

  #initialize static variables to pass to workers
  stabsel_init_param <- stabsel_init(listX = data_split_by_condition, nreps = nreps)


  #########################################################
  #**Initialize workers and perform stability selection **#
  #########################################################

  #print message to user
  message(paste0('Using Lambda hyper-parameter: ', optimized_lambda,'!'))
  message(paste0("stabilitySelection will be performed with ", nreps, " replicates"))

  if(subSample){

    message("Additional sub-sampling will be performed on uneven groups")
    ss_function <- "CGM_AHP_stabsel_subsample"
  }else if(!subSample){

    message("No additional sub-sampling will be performed. Sample groups will both be randomly
            sampled 50%")
    ss_function <- "CGM_AHP_stabsel"
  }

  stab_sel <- BiocParallel:: bplapply(X = 1:nreps,
                                      FUN = ss_function,
                                      init_param = stabsel_init_param,
                                      listX = data_split_by_condition,
                                      lastar = optimized_lambda,
                                      BPPARAM = BPPARAM,
                                      BPOPTIONS = bpoptions(progressbar = TRUE, tasks = 10))


  # #set progress bar
  # pbapply::pboptions(type = "timer", char = c('='), style = 5)
  #
  # #print lambda used
  # message(paste0('Using Lambda hyper-parameter: ', optimized_lambda,'!'))
  #
  # if(runParallel){
  #
  #   message('stabilitySelection() will be run in parallel ...', appendLF = TRUE)
  #
  #   #create independent processes to run reps in parallel
  #   cl <- parallel::makeCluster(nCores, type = 'PSOCK')
  #   on.exit(stopCluster(cl))
  #
  #   #pass necessary objects to independent workers
  #   parallel::clusterExport(cl = cl, varlist = c("CGM_AHP_tune",
  #                                                "CGM_AHP_train",
  #                                                "CGM_AHP_stabsel_subsample",
  #                                                "CGM_AHP_stabsel",
  #                                                "matTr"))
  #   parallel::clusterEvalQ(cl = cl, c(library("MASS"), library("glasso"), library("Matrix")))
  #
  #   if (subSample){
  #     message("Stability selection WITH additional subsampling using Guo et al ...\n", appendLF = TRUE)
  #
  #     stab_sel <- pbapply::pblapply(cl = cl,
  #                                   FUN = "CGM_AHP_stabsel_subsample",
  #                                   X = 1:nreps,
  #                                   init_param = stabsel_init_param,
  #                                   listX = data_split_by_condition,
  #                                   lastar = optimized_lambda)
  #   } else {
  #     message("Stability selection WITHOUT additional subsampling using Guo et al ...\n", appendLF = TRUE)
  #
  #     stab_sel <- pbapply::pblapply(cl = cl,
  #                                   FUN ="CGM_AHP_stabsel",
  #                                   X = 1:nreps,
  #                                   init_param = stabsel_init_param,
  #                                   listX = data_split_by_condition,
  #                                   lastar = optimized_lambda)
  #   }
  # }else{
  #
  #   message('stabilitySelection() will be run sequentially ...', appendLF = TRUE)
  #
  #   if (subSample){
  #     message("Stability selection WITH additional subsampling using Guo et al ...\n", appendLF = TRUE)
  #
  #     stab_sel <- pbapply::pblapply(FUN = "CGM_AHP_stabsel_subsample",
  #                                   X = 1:nreps,
  #                                   init_param = stabsel_init_param,
  #                                   listX = data_split_by_condition,
  #                                   lastar = optimized_lambda)
  #   } else {
  #     message("Stability selection WITHOUT additional subsampling using Guo et al ...\n", appendLF = TRUE)
  #
  #     stab_sel <- pbapply::pblapply(FUN ="CGM_AHP_stabsel",
  #                                   X = 1:nreps,
  #                                   init_param = stabsel_init_param,
  #                                   listX = data_split_by_condition,
  #                                   lastar = optimized_lambda)
  #   }
  # }

  #add empty line after progress bar
  message("", appendLF = TRUE)

  #####################################
  #**Concatenate results for output **#
  #####################################

  #initiate list for stability selection raw results
  selection_results <- vector("list", length(networkGroups(object)))
  names(selection_results) <- networkGroups(object)

  #initiate list for stability selection results converted to probabilities
  selection_probabilities <- vector("list", length(networkGroups(object)))
  names(selection_probabilities) <- networkGroups(object)


  #reduce results to one matrix and calculate selection probabilities
  for (k in 1:length(selection_results)){
    selection_results[[k]] <- lapply(stab_sel, function(r) r$mat[[k]])
    selection_results[[k]] <- Reduce("+", selection_results[[k]])

    if (subSample){
      message(paste0("Calculating selection probabilities WITH subsampling for...", names(selection_results)[[k]],"..."), appendLF = TRUE)
      selection_probabilities[[k]] <- selection_results[[k]]/(nreps)
    } else {
      message(paste0("Calculating selection probabilities WITHOUT subsampling for...",names(selection_results)[[k]],"..."), appendLF = TRUE)
      selection_probabilities[[k]] <- selection_results[[k]]/(2 * nreps)
    }
  }

  selectionResults(object) <- selection_results
  selectionProbabilities(object) <- selection_probabilities

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
#' @param optimal_lambda Lambda hyperparameter to be used in analysis. Not necessary if BICtune()
#'        or stabilitySelection() were already run.
#' @param eps A significance cut-off for thresholding the adjacency matrix. This can be left as default
#'
#' @return A DNEA object containing an adjacency matrix for the data network,
#'          determined by the glasso model
#'
#' @import zoo
#' @import glasso
#' @import glmnet
#' @import corpcor
#' @importFrom gdata lowerTriangle
#' @importFrom utils combn
#' @include JSEM.R
#' @include utilities.R
#' @export
getNeworks <- function(object,
                       optimal_lambda,
                       eps_threshold = 1e-06){

  ############################################
  #**Prepare data and initialize parameters**#
  ############################################

  num_samples <- numSamples(object)
  num_features <- numFeatures(object)

  Ip <- diag(rep(1,num_features))

  ##set up output to save the weighted and unweighted adjacency matrices from each model
  #weighted
  weighted_adjacency_matrices <- vector("list", length(networkGroups(object)))
  names(weighted_adjacency_matrices) <- networkGroups(object)

  #unweighted
  unweighted_adjacency_matrices <- vector("list", length(weighted_adjacency_matrices))
  names(unweighted_adjacency_matrices) <- names(weighted_adjacency_matrices)

  #separate the data by condition
  data_split_by_condition <- split_by_condition(dat = expressionData(object, type = "normalized"),
                                                condition_levels = networkGroups(object),
                                                condition_by_sample = networkGroupIDs(object))

  # stabilitySelection requires lambda hyper-parameter. Will use optimal_lambda if
  # supplied, otherwise looks for @hyperparameter[["optimized_lambda"]] in DNEAobject
  #
  # choosing lambda follows the following algorithm:
  # 1. if @hyperparamater[['optimized_lambda']] and optimal_lambda provided, optimal_lambda is used
  #    for analysis
  # 2. if @hyperparamater[['optimized_lambda']] and optimal_lambda missing, use
  #    @hyperparamater[['optimized_lambda']]
  # 3. if @hyperparamater[['optimized_lambda']] and optimal_lambda provided, optimal_lambda is used
  #    for analysis
  # 4. if both @hyperparamater[['optimized_lambda']] and optimal_lambda are missing, throw error.
  if(!is.null(optimizedLambda(object))){
    if(missing(optimal_lambda) == FALSE){

      optimized_lambda <- optimal_lambda
      warning('optimal_lambda argument was provided even though @hyperparameter[["optimized_lambda"]]
            already exists - optimal_lambda will be used in analysis')

    } else{

      optimized_lambda <- optimizedLambda(object)
    }
  } else{
    if(missing(optimal_lambda) == FALSE){

      optimized_lambda <- optimal_lambda
      optimizedLambda(object) <- optimal_lambda
      message('@hyperparameter[["optimized_lambda"]] was previously empty and now set to optimal_lambda argument')

    } else{

      # setting optimized_lambda = NULL will default to a lambda of sqrt(log(# features) / # samples)
      # in adjDGlasso_minimal
      optimized_lambda = NULL

      stop('No lambda value was supplied for the model - sqrt(log(# features) / # samples) will be
      used in the analyis. However, We highly recommend optimizing the lambda parameter by running
      BICtune(), or providing a calibrated lambda value using the optimal_lambda parameter prior to
           analysis.')

    }
  }

  #print lambda used
  message(paste0('Using Lambda hyper-parameter: ', optimized_lambda,'!'))

  #model will used selection weights based on stability selection if provided
  if (!is.null(selectionProbabilities(object))){

    model_weight_values <- lapply(selectionProbabilities(object),
                                  function(x) as.matrix(1/(1e-04 + x)))

    message('selection_probabilites from stability selection will be used in glasso model!')

  } else{

    message('No selection_probabilities were found. We recommend running
            stabilitySelection() prior to estimating the glasso model!')

    model_weight_values <- list(matrix(rep(1, num_features^2), num_features, num_features),
                                matrix(rep(1, num_features^2), num_features, num_features))

  }

  #add names to model weights list
  names(model_weight_values) <- names(selectionProbabilities(object))

  #############################################
  #**Estimate the partial correlation matrix**#
  #############################################

  for (k in networkGroups(object)){

    message(paste0('Estimating model for ', k, ' ...'), appendLF = TRUE)

    #fit the networks
    fit <- adjDGlasso_minimal(t(data_split_by_condition[[k]]),
                              weights= model_weight_values[[k]],
                              lambda = optimized_lambda)

    #grab the adjacency matrices
    weighted_adjacency_matrices[[k]] <- matrix(data = fit$Theta.glasso,
                                               nrow = 144, ncol = 144,
                                               dimnames = list(featureNames(object),
                                                               featureNames(object)))


  }

  ## threshold the weighted_adjacency_matrix as specified

  ## Get the unweighted adjacency matrix by thresholding the partial correlations
  for (k in names(weighted_adjacency_matrices)){
    unweighted_adjacency_matrices[[k]] <- abs(weighted_adjacency_matrices[[k]]) >= matrix(rep(eps_threshold, num_features^2), num_features, num_features)
  }

  #####################################
  #**Concatenate results for output **#
  #####################################

  #print message for total edges
  message(paste0("Number of edges in ", names(unweighted_adjacency_matrices)[[1]],": ", sum(unweighted_adjacency_matrices[[1]])/2), appendLF = TRUE)
  message(paste0("Number of edges in ", names(unweighted_adjacency_matrices)[[2]],": ", sum(unweighted_adjacency_matrices[[2]])/2), appendLF = TRUE)

  #store the adjacency matrices in DNEAresults object
  adjacencyMatrix(x = object, weighted = TRUE) <- weighted_adjacency_matrices
  adjacencyMatrix(x = object, weighted = FALSE) <- unweighted_adjacency_matrices

  #######################
  #**Create Edge List **#
  #######################

  #initiate output dataframe
  pairs <- combn(as.character(featureNames(object)), 2, simplify=FALSE)
  edge_list <- data.frame(Metabolite.A=rep(0,length(pairs)), Metabolite.B=rep(0,length(pairs)),
                          pcor.0=rep(0,length(pairs)), pcor.1=rep(0,length(pairs)),
                          check.names = FALSE)

  #concatenate results into dataframe
  edge_list[,1:2] <- do.call(rbind, pairs)
  edge_list[,3] <- lowerTriangle(weighted_adjacency_matrices[[1]])
  edge_list[,4] <- lowerTriangle(weighted_adjacency_matrices[[2]])
  edge_list$edge <- rep(NA, length(pairs)) #non-edge
  edge_list$edge[which((abs(edge_list$pcor.0) >= eps_threshold)*(abs(edge_list$pcor.1) >= eps_threshold) == 1)] <- "Both" #common edge
  edge_list$edge[which((abs(edge_list$pcor.0) >= eps_threshold)*(abs(edge_list$pcor.1)  < eps_threshold) == 1)] <- names(weighted_adjacency_matrices)[[1]]
  edge_list$edge[which((abs(edge_list$pcor.0) <  eps_threshold)*(abs(edge_list$pcor.1) >= eps_threshold) == 1)] <- names(weighted_adjacency_matrices)[[2]]
  edge_list <- edge_list[!is.na(edge_list$edge),]
  rownames(edge_list) <- NULL

  edgeList(object) <- edge_list
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
#'
#' @return A DNEAobject containing sub-network determinations for the nodes within the input network.
#'        A summary of the consensus clustering results can be viewed using getClusterResults().
#'        Sub-network classification for each node can be found in the node_list slot of the returned
#'        DNEAobject.
#'
#' @import igraph
#' @include preprocess_lib.R
#' @include utilities.R
#' @export
runConsensusCluster <- function(object, tau = 0.5){


  ###########################
  #**tau must be above 0.5**#
  ###########################

  if(tau < 0.5 | tau > 1.0) stop(paste0('tau corresponds to a percent agreement among the clustering methods. ',
                                        ' As such, tau must be greater than 0.5 and less than 1! ',
                                        'Clustering results below this threshold are not reliable -',
                                        'Please see user documentation for more information!'))

  #####################################
  #**Join the two condition networks**#
  #####################################

  #create list to hold graph from adjacency matrix
  adjacency_matrix_graphs <- vector("list", length(adjacencyMatrix(object, weighted = TRUE)))
  names(adjacency_matrix_graphs) <- names(adjacencyMatrix(object, weighted = TRUE))

  for (loop_el in names(adjacencyMatrix(object, weighted = TRUE))) {

    adjacency_graph <- graph_from_adjacency_matrix(adjacencyMatrix(object, weighted = TRUE)[[loop_el]], mode="undirected", weighted = TRUE)
    V(adjacency_graph)$name <- as.character(featureNames(object))
    adjacency_matrix_graphs[[loop_el]] <- adjacency_graph
  }

  #join adjacency matrix graphs
  joint_graph <- igraph::union(adjacency_matrix_graphs[[1]], adjacency_matrix_graphs[[2]])

  #modify the joint_graph
  jointLayout <- layout_nicely(joint_graph)
  E(joint_graph)$lty <- 1
  E(joint_graph)$color <- "black"
  E(joint_graph)$lty[is.na(E(joint_graph)$weight_2)] <- 2 #Group_1
  E(joint_graph)$lty[is.na(E(joint_graph)$weight_1)] <- 3 #Group_2
  E(joint_graph)$color[is.na(E(joint_graph)$weight_2)] <- "green" #Group_1
  E(joint_graph)$color[is.na(E(joint_graph)$weight_1)] <- "red"   #Group_2

  if(!is.null(nodeList(object)[, "DEstatus"])){
    V(joint_graph)$color <- ifelse(nodeList(object)[, "DEstatus"], "purple", "white")
    V(joint_graph)$DE <- nodeList(object)[, "DEstatus"]
  } else{
    V(joint_graph)$DE <- rep(NA, numFeatures(object))
  }

  ###########################################################
  #**ensemble community detection with consusus clustering**#
  ###########################################################

  #run consensus cluster algorithm
  # fit <- run_consensus_cluster(joint_graph,tau=tau0,method="ensemble", runParallel = runParallel, nCores = nCores)
  fit <- run_consensus_cluster(joint_graph, tau=tau)
  consensus_membership <- fit$final_consensus_cluster

  #initiate output matrix
  subnetwork_results <- matrix(0, nrow=length(unique(consensus_membership)), numFeatures(object),
                               dimnames = list(paste0("subnetwork",1:length(unique(consensus_membership))),
                                               sapply(1:length(joint_graph), function(x) names(joint_graph[[x]]))))

  #gather results
  for (j in 1:nrow(subnetwork_results)){
    subnetwork_results[j, which(consensus_membership == j)] <- 1
  }


  #####################################
  #**Concatenate results for output **#
  #####################################

  summary_list <- list()
  for (loop_cluster in 1:nrow(subnetwork_results) ){
    cluster_c <- induced.subgraph(joint_graph, V(joint_graph)$name[(subnetwork_results[loop_cluster,] == 1)])
    summary_list[[loop_cluster]] <- data.frame("number_of_nodes"=length(V(cluster_c)),
                                               "number_of_edges"=length(E(cluster_c)),
                                               "number_of_DE.nodes"=sum(as.numeric(table(V(cluster_c)$DE)[names(table(V(cluster_c)$DE))==TRUE])),
                                               "number_of_DE.edges"=sum(as.numeric(table(E(cluster_c)$color)[names(table(E(cluster_c)$color)) %in% c("red", "green")])),
                                               check.names = FALSE)
  }

  summary_stat <- data.frame("subnetworks"= rownames(subnetwork_results), do.call(rbind, summary_list), check.names = FALSE)

  nodeList(object)[["membership"]] <- consensus_membership
  object@consensus_clustering <- new(Class = "consensusClusteringResults",
                                     summary = summary_stat,
                                     subnetwork_membership = data.frame(subnetwork_results),
                                     adjacency_graphs = append(adjacency_matrix_graphs, list(joint_graph = joint_graph)))

  return(object)

}
#' Runs the NetGSA algorithm described in Hellstern et al.
#'
#' runNetGSA takes as input a DNEA object and utilizes the adjacency matrix calculated in getNetworks()
#' to perform NetGSA on the two networks to determine differences in subnetwork (pathway) expression.
#'
#' @param object A DNEAobject
#' @param min_size The minimum size of metabolic modules for enrichment analysis
#' @returns A DNEAobject containing containing results from NetGSA. Pathway expression differences for
#'          each node can be found in the node_list. A summary of the NetGSA results can be viewed
#'          using getNetGSAresults().
#'
#' @importFrom stats p.adjust
#' @import igraph
#' @import corpcor
#' @importFrom netgsa NetGSA
#' @include utilities.R
#' @export
runNetGSA <- function(object, min_size = 5){

  #################################
  #**Prepare data and run netGSA**#
  #################################

  ##set input variables
  adjacency_matrices <- list(list(adjacencyMatrix(x = object, weighted = TRUE)[[1]]),
                             list(adjacencyMatrix(x = object, weighted = TRUE)[[2]]))
  expression_data <- t(expressionData(object, type = "input"))
  data_groups <- ifelse(networkGroupIDs(object) == networkGroups(object)[1], 1, 2)
  subnetworks <- as.matrix(subnetworkMembership(object))

  ##filter subnetworks to only include those greater than or equal to min_size
  filtered_subnetworks <- subnetworkMembership(object)
  filtered_subnetworks <- as.matrix(filtered_subnetworks[rowSums(filtered_subnetworks) >= min_size, ])

  ##run netgsa
  netgsa_results <- NetGSA(A = adjacency_matrices,
                           x = expression_data,
                           group = data_groups,
                           pathways = filtered_subnetworks,
                           lklMethod = "REML",
                           minsize = min_size)



  # #separate the data by condition
  # separated_conditions_data <- split_by_condition(dat = expressionData(object, type = "input"),
  #                                                 condition_levels = networkGroups(object),
  #                                                 condition_by_sample = networkGroupIDs(object))

  # out.netgsa <- NetGSA(adjacencyMatrix(x = object, weighted = TRUE),
  #                      x = cbind(separated_conditions_data[[1]], separated_conditions_data[[2]]),
  #                      y = c(rep(1, ncol(separated_conditions_data[[1]])), rep(2, ncol(separated_conditions_data[[2]]))),
  #                      B = as.matrix(filtered_subnetworks), lklMethod = "REML")
  #####################################
  #**Concatenate results for output **#
  #####################################

  #add netGSA results to Node list
  nodeList(object)[["mean1"]] <- as.vector(netgsa_results$beta[[1]])
  nodeList(object)[["mean2"]] <- as.vector(netgsa_results$beta[[2]])
  nodeList(object)[["meanchange"]] <- netgsa_results$beta[[2]] - netgsa_results$beta[[1]]
  nodeList(object)[["mc.notes"]] <- paste(networkGroups(object)[[2]], 'over', networkGroups(object)[[1]])

  #concatenate netGSA summary output
  res <- data.frame(CCsummary(object)[CCsummary(object)$number_of_nodes >= min_size, ],
                    "NetGSA_pval"= netgsa_results$results$pval,
                    "NetGSA_pFDR"= netgsa_results$results$pFdr,
                    check.names = FALSE)

  #order netGSA results by FDR
  res <- res[order(res$NetGSA_pFDR),]
  rownames(res) <- 1:nrow(res)

  #Change nodelist membership to be indicative of new order
  # nodeList(object)[["membership"]][!(nodeList(object)[["membership"]] %in% gsub('subnetwork','',res$subnetworks))] <- NA
  # nodeList(object)[["membership"]] <- rownames(res)[match(nodeList(object)[["membership"]], as.numeric(gsub('subnetwork','',res$subnetworks)))]
  # nodeList(object)[["membership"]] <- as.numeric(nodeList(object)[["membership"]])

  #change res rownames to match new order
  res$subnetworks <- paste0("subnetwork ",rownames(res))

  #update DNEAobject
  netGSAresults(object) <- res

  return(object)
}













