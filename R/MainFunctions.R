#' @include JSEM.R
#' @include utilities.R
#' @include preprocess_lib.R
#' @include all-methods.R
#' @include all-generics.R
#' @include all-classes.R
#'
NULL

#' Optimize the lambda regularization parameter for the glasso-based network models using Bayesian-information Criterion
#'
#' This function will calculate the Bayesian information criterion (BIC) and liklihood for a range of lambda values
#' determined by the number of features within the dataset. The lambda value with the minimum BIC score is
#' the optimal lambda value for the respective dataset and is stored in the DNEAresults object for
#'  use in stability selection using \code{\link{stabilitySelection}} and network generation using
#'  \code{\link{getNetworks}}
#'
#' @param object A DNEAresults object. See \code{\link{createDNEAobject}}
#' @param lambda_values **OPTIONAL** A list of values to test while optimizing the lambda parameter.
#'  If not provided, a set of lambda values are chosen based on the theoretical value for the
#'   asymptoticly valid lambda. More information about this can be found in the details section
#' @param eps_threshold A significance cut-off for thresholding network edges
#'        The default value is 1e-06. This value generally should not change.
#' @param eta_value default parameter ??. Default is 0.1
#' @param BPPARAM A BiocParallel object
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{createDNEAobject}}
#'
#' @references Guo J, Levina E, Michailidis G, Zhu J. Joint estimation of multiple graphical models. Biometrika. 2011 Mar;98(1):1-15. doi: 10.1093/biomet/asq060. Epub 2011 Feb 9. PMID: 23049124; PMCID: PMC3412604.
#'
#' @return A DNEAresults object containing the BIC and liklihood scores for every lambda value tested, as well as
#'         the optimized lambda value
#'
#' @examples
#' #import example data
#' data(TEDDY)
#'
#' #initiate DNEAresults object
#' DNEA <- createDNEAobject(expression_data = TEDDY,
#'                          project_name = "TEDDYmetabolomics",
#'                          case = "DM:case",
#'                          control = "DM:control")
#'
#' #optimize lambda parameter
#' DNEA <- BICtune(object = DNEA, BPPARAM = bpparam())
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

  #check valid object
  validObject(object)

  return(object)
}

#' Stability selection to calculate selection probabilities for every possible feature-feature interaction within the data
#'
#' This function randomly samples 50% of samples for each condition a number of times specified by nreps
#' and fits a glasso model with the sampled data. A feature-feature interaction is considered present if the
#' partial correlation value is above 1e-5. The resulting adjacency matrices are summed together and
#' selection probabilities for each feature-feature interaction are calculated by dividing the number of replicates
#' performed by 2 (two models are fit with each half of the randomly sampled data, respectively, per replicate).
#'
#' The subSample parameter should be set to TRUE if the sample numbers across condition are uneven - this will
#' ensure that samples are randomly sampled in a way that ensures equal representation in stability selection. The
#' principles of stability selection remain similar with both methods, however, there are a few caveats. When
#' subSample = TRUE is set, The sample groups are stabilized by randomly sampling a smaller portion of the bigger group.
#' Due to the smaller sample sizes per sampling, only one model is fit per replicate.
#'
#' @param object A DNEAresults object
#' @param subSample A boolean that specifies whether the number of samples are unevenly split
#'         by condition and, therefore, should be adjusted for when randomly sampling.
#' @param nreps The total number of replicates to perform stability selection; The default is 500.
#'        Performing a high number of replicates becomes increasingly important as datasets decrease in size.
#' @param optimal_lambda ***OPTIONAL*** The optimal lambda value to be used in the model. This parameter is only
#'        necessary if BICtune() is not performed
#' @param BPPARAM a BiocParallel object
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{createDNEAobject}}, \code{\link{BICtune}}
#'
#' @references Ma J, Karnovsky A, Afshinnia F, Wigginton J, Rader DJ, Natarajan L, Sharma K, Porter AC, Rahman M, He J, Hamm L, Shafi T, Gipson D, Gadegbeku C, Feldman H, Michailidis G, Pennathur S. Differential network enrichment analysis reveals novel lipid pathways in chronic kidney disease. Bioinformatics. 2019 Sep 15;35(18):3441-3452. doi: 10.1093/bioinformatics/btz114. PMID: 30887029; PMCID: PMC6748777. \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6748777/}
#'
#' @return A DNEAresults object after populating the stable_networks slot of the object. It contains the selection
#' results from stability selection as well as the calculated selection probabilities
#'
#' @examples
#' #import example data
#' data(TEDDY)
#'
#' #initiate DNEAresults object
#' DNEA <- createDNEAobject(expression_data = TEDDY,
#'                          project_name = "TEDDYmetabolomics",
#'                          case = "DM:case",
#'                          control = "DM:control")
#'
#' #optimize lambda parameter
#' DNEA <- BICtune(object = DNEA, BPPARAM = bpparam())
#'
#' # perform stability selection
#' DNEA <- stabilitySelection(object = DNEA, subSample = FALSE, nreps = 5, BPPARAM = bpparam())
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
                               nreps = 500,
                               optimal_lambda,
                               BPPARAM = bpparam()){

  # stabilitySelection requires lambda hyper-parameter. Will use optimal_lambda if
  # supplied, otherwise looks for @hyperparameter[["optimized_lambda"]] in DNEAobject
  #
  # choosing lambda follows the following algorithm:
  # 1. if @hyperparamater[['optimized_lambda']] and optimal_lambda provided, optimal_lambda is used
  #    for analysis
  # 2. if optimal_lambda provided and @hyperparamater[['optimized_lambda']] missing, use optimal lambda
  #    and set to store as new lambda value
  # 3. if @hyperparamater[['optimized_lambda']] provided and optimal_lambda missing, @hyperparamater[['optimized_lambda']]
  #    is used for analysis
  # 4. if both @hyperparamater[['optimized_lambda']] and optimal_lambda are missing, default to sqrt(log(p)/n) and give warning
  if(!missing(optimal_lambda)){
    if(!is.null(optimizedLambda(object))){

      optimized_lambda <- optimal_lambda
      warning('optimal_lambda argument was provided even though @hyperparameter[["optimized_lambda"]]
            already exists - optimal_lambda will be used in analysis')
    }else{

      optimized_lambda <- optimal_lambda
      optimizedLambda(object) <- optimal_lambda
      message('@hyperparameter[["optimized_lambda"]] was previously empty and now set to optimal_lambda argument')
    }
  }else if(!is.null(optimizedLambda(object))){

    optimized_lambda <- optimizedLambda(object)
  }else{

    # setting optimized_lambda = NULL will default to a lambda of sqrt(log(# features) / # samples)
    # in adjDGlasso_minimal
    optimized_lambda = NULL

    stop('No lambda value was supplied for the model - sqrt(log(# features) / # samples) will be
      used in the analyis. However, We highly recommend optimizing the lambda parameter by running
      BICtune(), or providing a calibrated lambda value using the optimal_lambda parameter prior to
           analysis.')
  }

  #split data by condition
  data_split_by_condition = lapply(split_by_condition(dat = expressionData(object, type = "normalized"),
                                    condition_levels = networkGroups(object),
                                    condition_by_sample = networkGroupIDs(object)), function(d) t(d))

  #initialize static variables to pass to workers
  stabsel_init_param <- stabsel_init(listX = data_split_by_condition, nreps = nreps)


  #check that groups are sufficiently uneven if subSample selected
  if(subSample){
    if((1.3* stabsel_init_param[["min_num_sample"]]) > max(stabsel_init_param[["num_samples"]])) stop(paste0("The condition groups are not sufficiently uneven to randomly sample apropriately.\n",
                                                                                                             "Please perform stability selection WITHOUT additional sub-sampling"))
  }

  ##perform stability selection
  #print message to user
  message(paste0('Using Lambda hyper-parameter: ', optimized_lambda,'!'))
  message(paste0("stabilitySelection will be performed with ", nreps, " replicates"))

  if(subSample){

    #with additional sub-sampling
    message("Additional sub-sampling will be performed on uneven groups")
    ss_function <- "CGM_AHP_stabsel_subsample"
  }else if(!subSample){

    #without additional sub-sampling
    message("No additional sub-sampling will be performed. Sample groups will both be randomly sampled 50%")
    ss_function <- "CGM_AHP_stabsel"
  }

  #run SS
  stab_sel <- BiocParallel:: bplapply(X = seq(1, nreps),
                                      FUN = ss_function,
                                      init_param = stabsel_init_param,
                                      listX = data_split_by_condition,
                                      lastar = optimized_lambda,
                                      BPPARAM = BPPARAM,
                                      BPOPTIONS = bpoptions(progressbar = TRUE, tasks = 10))

  #add empty line after progress bar
  message("", appendLF = TRUE)

  ##concatenate results for output
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

  #check valid object
  validObject(object)

  return(object)
}

#' Construct the GLASSO-based biological Networks
#'
#' This function utilizes the lambda parameter tuned using \code{\link{BICtune}} *(may also be user-specified)*
#' and *optionally* the feature-feature edge selection probabilities using \code{\link{stabilitySelection}} to jointly estimate the biological networks
#' within the dataset.
#'
#' @param object A DNEAresults object
#' @param optimal_lambda ***OPTIONAL*** The lambda hyperparameter to be used in analysis. Not necessary if \code{\link{BICtune}}
#'        or \code{\link{stabilitySelection}} were already performed
#' @param eps_threshold A numeric value between 0 and 1 by which to threshold the partial correlation values for edge identification.
#' Edges with a partial correlation value below this threshold will be zero'd out from the adjacency matrix.
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{createDNEAobject}}, \code{\link{BICtune}}, \code{\link{stabilitySelection}},
#' \code{\link{edgeList}}, \code{\link{adjacencyMatrix}}
#'
#' @references
#' Ma J, Karnovsky A, Afshinnia F, Wigginton J, Rader DJ, Natarajan L, Sharma K, Porter AC, Rahman M, He J, Hamm L, Shafi T, Gipson D, Gadegbeku C, Feldman H, Michailidis G, Pennathur S. Differential network enrichment analysis reveals novel lipid pathways in chronic kidney disease. Bioinformatics. 2019 Sep 15;35(18):3441-3452. doi: 10.1093/bioinformatics/btz114. PMID: 30887029; PMCID: PMC6748777. \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6748777/}
#'
#' Iyer GR, Wigginton J, Duren W, LaBarre JL, Brandenburg M, Burant C, Michailidis G, Karnovsky A. Application of Differential Network Enrichment Analysis for Deciphering Metabolic Alterations. Metabolites. 2020; 10(12):479. \url{https://doi.org/10.3390/metabo10120479}
#'
#' @return A DNEAresults object after populating the adjaceny_matrix and edge_list slots with the corresponding
#' adjacency_matrix for each sample condition as well as the network edge list.
#'
#' @examples
#' #import example data
#' data(TEDDY)
#'
#' #initiate DNEAresults object
#' DNEA <- createDNEAobject(expression_data = TEDDY,
#'                          project_name = "TEDDYmetabolomics",
#'                          case = "DM:case",
#'                          control = "DM:control")
#'
#' #optimize lambda parameter
#' DNEA <- BICtune(object = DNEA, BPPARAM = bpparam())
#'
#' # perform stability selection
#' DNEA <- stabilitySelection(object = DNEA, subSample = FALSE, nreps = 5, BPPARAM = bpparam())
#'
#' #construct the networks
#' DNEA <- getNetworks(object = DNEA)
#'
#' #now we can plot the group networks
#' plotNetworks(object = DNEA, type = "group_networks")
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
getNetworks <- function(object,
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

  # getNetworks() requires lambda hyper-parameter. Will use optimal_lambda if
  # supplied, otherwise looks for @hyperparameter[["optimized_lambda"]] in DNEAobject
  #
  # choosing lambda follows the following algorithm:
  # 1. if @hyperparamater[['optimized_lambda']] and optimal_lambda provided, optimal_lambda is used
  #    for analysis
  # 2. if optimal_lambda provided and @hyperparamater[['optimized_lambda']] missing, use optimal lambda
  #    and set to store as new lambda value
  # 3. if @hyperparamater[['optimized_lambda']] provided and optimal_lambda missing, @hyperparamater[['optimized_lambda']]
  #    is used for analysis
  # 4. if both @hyperparamater[['optimized_lambda']] and optimal_lambda are missing, default to sqrt(log(p)/n) and give warning
  if(!missing(optimal_lambda)){
    if(!is.null(optimizedLambda(object))){

      optimized_lambda <- optimal_lambda
      warning('optimal_lambda argument was provided even though @hyperparameter[["optimized_lambda"]]
            already exists - optimal_lambda will be used in analysis')
    }else{

      optimized_lambda <- optimal_lambda
      optimizedLambda(object) <- optimal_lambda
      message('@hyperparameter[["optimized_lambda"]] was previously empty and now set to optimal_lambda argument')
    }
  }else if(!is.null(optimizedLambda(object))){

    optimized_lambda <- optimizedLambda(object)
  }else{

    # setting optimized_lambda = NULL will default to a lambda of sqrt(log(# features) / # samples)
    # in adjDGlasso_minimal
    optimized_lambda = NULL

    stop('No lambda value was supplied for the model - sqrt(log(# features) / # samples) will be
      used in the analyis. However, We highly recommend optimizing the lambda parameter by running
      BICtune(), or providing a calibrated lambda value using the optimal_lambda parameter prior to
           analysis.')
  }


  #print lambda used
  message(paste0('Using Lambda hyper-parameter: ', optimized_lambda,'!'))

  #model will used selection weights based on stability selection if provided
  if (!is.null(selectionProbabilities(object))){

    model_weight_values <- lapply(selectionProbabilities(object),
                                  function(x) as.matrix(1/(1e-04 + x)))

    message('selection_probabilites from stability selection will be used in glasso model!\n')

  } else{

    message("No selection_probabilities were found. We recommend running
            stabilitySelection() prior to estimating the glasso model!\n")

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

  ##input weighted_adjacency_matrices into object
  adjacencyMatrix(x = object, weighted = TRUE) <- weighted_adjacency_matrices

  ##filter the weighted_adjacency_matrices by eps_threshold and create unweighed_adjacency_matrices plus edge list
  object <- filterNetworks(data = object, pcor = eps_threshold)

  #check valid object
  validObject(object)

  return(object)
}
#' Identify metabolic modules within the biological networks using a consensus clustering approach
#'
#' This function clusters the biological networks constructed using \code{\link{getNetworks}} using a consensus clustering
#' approach described in Ma et al. *(Please see references for more details)*
#'
#' Seven clustering algorithms from the \code{\link{igraph}} package:
#' 1. \code{\link{igraph::cluster_edge_betweenness}}
#'
#' 2. \code{\link{igraph::cluster_fast_greedy}}
#'
#' 3. \code{\link{igraph::cluster_infomap}}
#'
#' 4. \code{\link{igraph::cluster_label_prop}}
#'
#' 5. \code{\link{igraph::cluster_louvain}}
#'
#' 6. \code{\link{igraph::cluster_walktrap}}
#'
#' 7. \code{\link{igraph::cluster_leading_eigen}}
#'
#' are performed iteratively on the adjacency matrix constructed using \code{\link{getNetworks}} until a consensus
#' is reached on resulting subnetwork membership, or the specified max_iterations is reached.
#'
#'
#' @param object A DNEAresults object
#' @param tau The % agreement threshold among the clustering algorithms for a node to be included in a subnetwork
#' @param max_iterations The maximum number of replicates of the clustering algorithms to perform before consensus is reached
#' @param eps_threshold A cut-off value thresholding the adjacency matrix. Edges with a partial correlation value below this
#' threshold will be removed.
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{createDNEAobject}}, \code{\link{BICtune}}, \code{\link{stabilitySelection}},
#' \code{\link{edgeList}}, \code{\link{adjacencyMatrix}}, \code{\link{getNetworks}}, \code{\link{edgeList}},
#' \code{\link{thresholdNetworks}}, \code{\link{plotNetworks}}
#'
#' @references
#' Ma J, Karnovsky A, Afshinnia F, Wigginton J, Rader DJ, Natarajan L, Sharma K, Porter AC, Rahman M, He J, Hamm L, Shafi T, Gipson D, Gadegbeku C, Feldman H, Michailidis G, Pennathur S. Differential network enrichment analysis reveals novel lipid pathways in chronic kidney disease. Bioinformatics. 2019 Sep 15;35(18):3441-3452. doi: 10.1093/bioinformatics/btz114. PMID: 30887029; PMCID: PMC6748777. \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6748777/}
#'
#' @return A DNEAobject containing sub-network determinations for the nodes within the input network.
#'        A summary of the consensus clustering results can be viewed using getClusterResults().
#'        Sub-network classification for each node can be found in the node_list slot of the returned
#'        DNEAobject.
#'
#' @examples
#' #import example data
#' data(TEDDY)
#'
#' #initiate DNEAresults object
#' DNEA <- createDNEAobject(expression_data = TEDDY,
#'                          project_name = "TEDDYmetabolomics",
#'                          case = "DM:case",
#'                          control = "DM:control")
#'
#' #optimize lambda parameter
#' DNEA <- BICtune(object = DNEA, BPPARAM = bpparam())
#'
#' # perform stability selection
#' DNEA <- stabilitySelection(object = DNEA, subSample = FALSE, nreps = 500, BPPARAM = bpparam())
#'
#' #construct the networks
#' DNEA <- getNetworks(object = DNEA)
#'
#' #view the edgelist
#' edgeList(object)
#'
#' #identify metabolic modules via consensus clustering
#' DNEA <- runConsensusCluster(object = DNEA)
#'
#' #we can also plot the subnetworks
#' plotNetworks(object = DNEA, type = "subnetworks", subnetwork = 1)
#'
#' @import igraph
#' @include preprocess_lib.R
#' @include utilities.R
#' @export
runConsensusCluster <- function(object,
                                tau = 0.5,
                                max_iterations = 5,
                                eps_threshold = 1e-06){


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
  fit <- run_consensus_cluster(joint_graph, tau=tau, max_iterations = max_iterations)
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

  #check valid object
  validObject(object)

  return(object)

}
#' Identify metabolic modules that are enriched across experimental conditions
#'
#' This function performs pathway enrichment analysis on the metabolic modules identified via \code{\link{runConsensusCluster}}
#' using the \code{\link{netgsa}} algorithm.
#'
#' @param object A DNEAobject
#' @param min_size The minimum size of metabolic modules for enrichment analysis
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{netGSAresults}}, \code{\link{plotNetworks}}, \code{\link{nodeList}}
#'
#' @references
#' Hellstern M, Ma J, Yue K, Shojaie A (2021) netgsa: Fast computation and interactive visualization for topology-based pathway enrichment analysis. PLoS Comput Biol 17(6): e1008979. \url{https://doi.org/10.1371/journal.pcbi.1008979}
#'
#'
#' @returns A DNEAobject containing containing results from NetGSA. Pathway expression differences for
#'          each node can be found in the node_list. A summary of the NetGSA results can be viewed
#'          using getNetGSAresults()
#'
#' @examples
#' #import example data
#' data(TEDDY)
#'
#' #initiate DNEAresults object
#' DNEA <- createDNEAobject(expression_data = TEDDY,
#'                          project_name = "TEDDYmetabolomics",
#'                          case = "DM:case",
#'                          control = "DM:control")
#'
#' #optimize lambda parameter
#' DNEA <- BICtune(object = DNEA, BPPARAM = bpparam())
#'
#' # perform stability selection
#' DNEA <- stabilitySelection(object = DNEA, subSample = FALSE, nreps = 5, BPPARAM = bpparam())
#'
#' #construct the networks
#' DNEA <- getNetworks(object = DNEA)
#'
#' #identify metabolic modules via consensus clustering
#' DNEA <- runConsensusCluster(object = DNEA)
#'
#' #perform pathway enrichment analysis using netGSA
#' DNEA <- runNetGSA(object = DNEA)
#'
#' #view the results
#' netGSAresults(DNEA)
#'
#' #save node and edge list for input to cytoscape
#' getNetworkFiles(DNEA)
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

  #update DNEAobject
  netGSAresults(object) <- res

  #check valid object
  validObject(object)

  return(object)
}













