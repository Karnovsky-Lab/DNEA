#' @include JSEM.R
#' @include utilities-external.R
#' @include utilities-internal.R
#' @include preprocess_lib.R
#' @include all-methods.R
#' @include all-generics.R
#' @include all-classes.R
#' @include start-here.R
NULL

#' Optimize the lambda regularization parameter for the glasso-based network models using Bayesian-information Criterion
#'
#' This function will calculate the Bayesian information criterion (BIC) and likelihood for a range of lambda values
#' that are automatically generated (\emph{please see \strong{Details} for more info}) or that are user-specified.
#' The lambda value with the minimum BIC score is the optimal lambda value for the dataset and is stored in the
#' DNEAresults object for use in stability selection using \code{\link{stabilitySelection}} and network generation using
#'  \code{\link{getNetworks}}
#'
#' @param object A \code{DNEAresults} object. See \code{\link{createDNEAobject}}
#' @param lambda_values **OPTIONAL** A list of values to test while optimizing the lambda parameter.
#'  If not provided, a set of lambda values are chosen based on the theoretical value for the
#'  asymptotically valid lambda. More information about this can be found in the details section
#' @param eps_threshold A significance cut-off for thresholding network edges.
#'        The default value is 1e-06. This value generally should not change.
#' @param eta_value default parameter ??. Default is 0.1
#' @param BPPARAM A \code{\link{BiocParallel}} object
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{createDNEAobject}}
#'
#' @references Guo J, Levina E, Michailidis G, Zhu J. Joint estimation of multiple graphical models. Biometrika. 2011 Mar;98(1):1-15. doi: 10.1093/biomet/asq060. Epub 2011 Feb 9. PMID: 23049124; PMCID: PMC3412604. \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3412604/}
#'
#' @details
#' There are several ways to optimize the lambda parameter for a glasso model - We utilize Bayesian-information criterion (BIC) to optimize the lambda parameter in
#' DNEAdev because it is a more balanced method and less computationally expensive. We can reduce the total number of values that
#' need to be tested in optimization by carefully selecting values around the asymptotically valid lambda for datasets with many samples and many features:
#'  \deqn{\lambda = \sqrt{ \ln (num. features) / num. samples}}{lambda = sqrt(ln(num. features) / num. samples)}
#' For smaller datasets, the asymptotically valid lambda is described by modifying the previous equation to include an unknown constant, c,
#' that needs to be determined mathematically. Therefore, to optimize lambda we modify the previous equation as follows:
#' \deqn{\lambda = c \sqrt{ \ln (num. features) / num. samples}}{lambda = c*sqrt(ln(num. features) / num. samples)}
#' where c takes on 15 evenly spaced values between 0.01 and 0.3. More information regarding the optimization method deployed here can be found
#' in the Guo et al. (2011) paper referenced below.
#'
#' @returns A \code{DNEAresults} object containing the BIC and likelihood scores for every lambda value tested, as well as
#'         the optimized lambda value
#'
#' @examples
#' #import BiocParallel package
#' library(BiocParallel)
#'
#' #import completed example data
#' data(TEDDYresults)
#'
#' #optimize lambda parameter
#' TEDDYresults <- BICtune(object = TEDDYresults, BPPARAM = bpparam())
#'
#' @import glasso
#' @importFrom BiocParallel bplapply bpparam bpoptions bptasks
#' @export
BICtune <- function(object,
                    lambda_values,
                    eps_threshold = 1e-06,
                    eta_value = 0.1,
                    BPPARAM = bpparam()){


  ##test for proper input
  if(!inherits(object, "DNEAresults")) stop('the input object should be of class "DNEAresults"!')

  ##prepare data
  dat <- split_by_condition(dat = expressionData(object, normalized = TRUE),
                            condition_levels = networkGroups(object),
                            condition_by_sample = networkGroupIDs(object))

  ##initialize input parameters
  n4cov <- max(vapply(dat, ncol, numeric(1)))
  trainX <- t(do.call(cbind, dat))
  trainY <- c(rep(1, ncol(dat[[1]])),
              rep(2, ncol(dat[[2]])))


  ##Pre-define a range of lambda values to evaluate during optimization if none are provided
  if(missing(lambda_values)){

    lambda_values <- seq(0.01, 0.3, 0.02)*sqrt(log(numFeatures(object))/n4cov)
  }else{

    lambda_values <- unlist(lambda_values)

    #check that lambda's are valid
    if(any(lambda_values < 0 | lambda_values > 1)) stop("The lambda parameter should be a value between 0 and 1 only!")
  }

  ##call internal tuning function to optimize lambda
  message("Optimizing the lambda hyperparameter using Bayesian-Information Criterion outlined in Guo et al. (2011)", appendLF = TRUE)
  message("A Link to this reference can be found in the function documentation by running ?BICtune() in the console", appendLF = TRUE)

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
  BIC_scores <- unlist(vapply(BIC_guo, function(a) a$BIC, numeric(1)))

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

  message("The optimal Lambda hyper-parameter has been set to: ", appendLF = FALSE)
  message(lastar_guo, appendLF = FALSE); message("!", appendLF = TRUE)

  #check valid object
  validObject(object)

  return(object)
}

#' Stability selection to calculate selection probabilities for every possible feature-feature interaction within the data
#'
#' This function randomly samples the input data and fits a glasso model with the sampled data for \strong{\emph{nreps}} number of replicates.
#' The resulting adjacency matrices are summed together and selection probabilities for each feature-feature interaction are calculated.
#' Stability selection is particularly useful for smaller datasets and when a large number of replicates are performed (the default is 500). The
#' exact method deployed varies slightly whether or not additional sub-sampling of the data is performed. More information can be
#' found in the \strong{\emph{Details}} section.
#'
#' @param object A \code{DNEAresults} object
#' @param subSample A boolean that specifies whether the number of samples are unevenly split
#'         by condition and, therefore, should be adjusted for when randomly sampling.
#' @param nreps The total number of replicates to perform in stability selection. The default is 500.
#' @param optimal_lambda \emph{OPTIONAL} - The optimal lambda value to be used in the model. This parameter is only
#'        necessary if \code{\link{BICtune}} is not performed
#' @param BPPARAM a BiocParallel object
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{createDNEAobject}}, \code{\link{BICtune}}
#'
#' @references Ma J, Karnovsky A, Afshinnia F, Wigginton J, Rader DJ, Natarajan L, Sharma K, Porter AC, Rahman M, He J, Hamm L, Shafi T, Gipson D, Gadegbeku C, Feldman H, Michailidis G, Pennathur S. Differential network enrichment analysis reveals novel lipid pathways in chronic kidney disease. Bioinformatics. 2019 Sep 15;35(18):3441-3452. doi: 10.1093/bioinformatics/btz114. PMID: 30887029; PMCID: PMC6748777. \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6748777/}
#'
#' Nicolai, M., & Peter, B. (2010). Stability selection. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 72(4), 417-473. \url{https://stat.ethz.ch/Manuscripts/buhlmann/stability.pdf}
#'
#' @details
#' Stability selection provides an additional approach by which to regularize the network model and create more robust results,
#' particularly when \strong{\emph{p >> n}}. Stability selection works by randomly sampling (without replacement) the input data many times and fitting
#' a glasso model to each subset of sampled data. The unwieghted adjacency matrix from each model is summed together
#' (A feature-feature interaction is considered present if the partial correlation value is above 1e-5), and the probability of
#' an edge being selected in a random subset of the data is calculated by dividing the number of times an edge was selected in
#' the replicates over the total number of replicates. This results in a selection probability for every possible feature-feature interaction
#' that is used to modify the regularization parameter via the following equation:
#' \deqn{\rho = \lambda*(1 / (0.0001 + selection.probability))}{ rho = lambda*(1 / (0.0001 + selection.probability))}
#'
#' However, when the sample groups are very unbalanced, randomly sampling strongly favors the larger group, resulting in
#' over representation of the aforementioned group. In order to combat this, setting subSample = TRUE modifies the random
#' sample by sub-sampling the groups individually to even out the numbers. In this method, 90% of the smaller group is randomly
#' sampled without replacement, and an additional 10% is randomly sampled without replacement from the entire group to preserve the
#' variance. The larger group is randomly sampled to have 1.3 times the number of samples present in the smaller group. This method ensures that
#' each group is equally represented in stability selection.\cr
#'
#' The principles of stability selection remain similar with both methods, however, there are a few small differences.
#' Stability selection \emph{without} additional sub-sampling randomly samples 50% of each group (without replacement) and fits a
#' model for both halves of the sampled data. Since nearly all of the data for the smaller group is used \emph{with} additional sub-sampling, only one
#' model is fit per replicate when subSample = TRUE. This means that at the default value of nreps = 500, 1000 randomly sampled
#' models are fit in total \emph{without} sub-sampling, but 500 randomly sampled models are fit in total \emph{with} sub-sampling.
#' More details about the stability approach deployed in this function can be found in Ma et al. (2019) referenced below.
#'
#'
#' @returns A \code{DNEAresults} object after populating the stable_networks slot of the object. It contains the selection
#' results from stability selection as well as the calculated selection probabilities.
#'
#' @examples
#' #import BiocParallel package
#' library(BiocParallel)
#'
#' #import completed example data
#' data(TEDDYresults)
#'
#' # perform stability selection
#' TEDDYresults <- stabilitySelection(object = TEDDYresults,
#'                                    subSample = FALSE,
#'                                    nreps = 4,
#'                                    BPPARAM = bpparam())
#'
#' @import glasso
#' @importFrom BiocParallel bplapply bpparam bpoptions bptasks
#' @export
stabilitySelection <- function(object,
                               subSample = FALSE,
                               nreps = 500,
                               optimal_lambda,
                               BPPARAM = bpparam()){

  ##test for proper input
  if(!inherits(object, "DNEAresults")) stop('the input object should be of class "DNEAresults"!')
  if(nreps < 1 | !is.numeric(nreps)) stop("nreps specifies the number of stability selection replicates to perform and should be an number greater than zero!")

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

    #check that lambda's are valid
    if(optimal_lambda < 0 | optimal_lambda > 1) stop("The lambda parameter should be a value between 0 and 1 only!")
    if(!is.null(optimizedLambda(object))){

      optimized_lambda <- optimal_lambda
      warning('optimal_lambda argument was provided even though @hyperparameter[["optimized_lambda"]] already exists",
              " - optimal_lambda will be used in analysis')
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
    optimized_lambda <- NULL

    warning("No lambda value was supplied for the model - sqrt(log(# features) / # samples) will beused in the analyis. ",
            "However, We highly recommend optimizing the lambda parameter by running BICtune(), ",
            "or providing a calibrated lambda value using the optimal_lambda parameter prior to analysis.")
  }

  #split data by condition
  data_split_by_condition <- lapply(split_by_condition(dat = expressionData(object, normalized = TRUE),
                                    condition_levels = networkGroups(object),
                                    condition_by_sample = networkGroupIDs(object)), function(d) t(d))

  #initialize static variables to pass to workers
  stabsel_init_param <- stabsel_init(listX = data_split_by_condition, nreps = nreps)


  #check that groups are sufficiently uneven if subSample selected
  if(subSample){
    if((1.3* stabsel_init_param[["min_num_samples"]]) > max(stabsel_init_param[["num_samples"]])) stop("The condition groups are not sufficiently uneven to randomly sample apropriately.\n",
                                                                                                       "Please perform stability selection WITHOUT additional sub-sampling")
  }

  ##perform stability selection
  #print message to user
  message("Using Lambda hyper-parameter: ", appendLF = FALSE)
  message(optimized_lambda, appendLF = FALSE); message("!", appendLF = TRUE)

  message("stabilitySelection will be performed with ", appendLF = FALSE)
  message(nreps, appendLF = FALSE);message(" replicates", appendLF = TRUE)

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
  for (k in seq(1, length(selection_results))){
    selection_results[[k]] <- lapply(stab_sel, function(r) r$mat[[k]])
    selection_results[[k]] <- Reduce("+", selection_results[[k]])

    if (subSample){

      message("Calculating selection probabilities WITH subsampling for...", appendLF = FALSE)
      message(names(selection_results)[[k]], appendLF = FALSE);message("...", appendLF = TRUE)

      selection_probabilities[[k]] <- selection_results[[k]]/(nreps)
    } else {

      message("Calculating selection probabilities WITHOUT subsampling for...", appendLF = FALSE)
      message(names(selection_results)[[k]], appendLF = FALSE);message("...", appendLF = TRUE)

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
#' This function constructs the biological network for each experimental condition using the joint estimation method described
#' in Ma et al. (2019) (\emph{please see references below}). If \code{\link{BICtune}} and \code{\link{stabilitySelection}}
#' were already run, the optimized lambda and selection probabilities from each function, respectively, will be used to
#' add regularization when constructing the networks (please see the \strong{\emph{Details}} section of \code{\link{stabilitySelection}}
#' for more information). Otherwise, \deqn{\lambda = \sqrt{\ln (num. features) / num. samples}}{ lambda = sqrt(ln(num. features) / num. samples)}
#' will be used as the regularization parameter.
#'
#' @param object A \code{DNEAresults} object
#' @param optimal_lambda \emph{OPTIONAL} - The lambda value to be used in analysis. Not necessary if \code{\link{BICtune}}
#'        or \code{\link{stabilitySelection}} were already performed
#' @param eps_threshold A numeric value between 0 and 1 by which to threshold the partial correlation values for edge identification.
#' Edges with an absolute partial correlation value below this threshold will be zero'd out from the adjacency matrix.
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{createDNEAobject}}, \code{\link{BICtune}}, \code{\link{stabilitySelection}},
#' \code{\link{edgeList}}, \code{\link{adjacencyMatrix}}
#'
#' @references
#' Ma J, Karnovsky A, Afshinnia F, Wigginton J, Rader DJ, Natarajan L, Sharma K, Porter AC, Rahman M, He J, Hamm L, Shafi T, Gipson D, Gadegbeku C, Feldman H, Michailidis G, Pennathur S. Differential network enrichment analysis reveals novel lipid pathways in chronic kidney disease. Bioinformatics. 2019 Sep 15;35(18):3441-3452. doi: 10.1093/bioinformatics/btz114. PMID: 30887029; PMCID: PMC6748777. \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6748777/}
#'
#' Iyer GR, Wigginton J, Duren W, LaBarre JL, Brandenburg M, Burant C, Michailidis G, Karnovsky A. Application of Differential Network Enrichment Analysis for Deciphering Metabolic Alterations. Metabolites. 2020 Nov 24;10(12):479. doi: 10.3390/metabo10120479. PMID: 33255384; PMCID: PMC7761243. \url{https://pubmed.ncbi.nlm.nih.gov/33255384/}
#'
#' @returns A \code{DNEAresults} object after populating the adjaceny_matrix and edge_list slots with the corresponding
#' adjacency_matrix for each sample condition as well as the network edge list.
#'
#' @examples
#' #import completed example data
#' data(TEDDYresults)
#'
#' #construct the networks
#' TEDDYresults <- getNetworks(object = TEDDYresults)
#'
#' #now we can plot the group networks
#' plotNetworks(object = TEDDYresults, type = "group_networks")
#'
#' @importFrom gdata lowerTriangle
#' @importFrom utils combn
#' @export
getNetworks <- function(object,
                       optimal_lambda,
                       eps_threshold = 1e-06){

  ##initialize input parameters
  num_samples <- numSamples(object)
  num_features <- numFeatures(object)
  Ip <- diag(rep(1, num_features))

  ##initiate output data structures
  #weighted adjacency matrices list
  weighted_adjacency_matrices <- vector("list", length(networkGroups(object)))
  names(weighted_adjacency_matrices) <- networkGroups(object)

  #unweighted adjacency matrices list
  unweighted_adjacency_matrices <- vector("list", length(weighted_adjacency_matrices))
  names(unweighted_adjacency_matrices) <- names(weighted_adjacency_matrices)

  ##test for proper input
  if(!inherits(object, "DNEAresults")) stop('the input object should be of class "DNEAresults"!')
  if(eps_threshold <=0 | eps_threshold >= 1) stop("The partial correlation threshold should be between 0 and 1 only!")

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

    #check that lambda's are valid
    if(optimal_lambda < 0 | optimal_lambda > 1) stop("The lambda parameter should be a value between 0 and 1 only!")
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
    optimized_lambda <- NULL

    warning("No lambda value was supplied for the model - sqrt(log(# features) / # samples) will beused in the analyis. ",
            "However, We highly recommend optimizing the lambda parameter by running BICtune(), ",
            "or providing a calibrated lambda value using the optimal_lambda parameter prior to analysis.")
  }
  #print lambda used
  message("Using Lambda hyper-parameter: ", appendLF = FALSE)
  message(optimized_lambda, appendLF = FALSE);message("!", appendLF = TRUE)

  ##separate the data by condition
  data_split_by_condition <- split_by_condition(dat = expressionData(object, normalized = TRUE),
                                                condition_levels = networkGroups(object),
                                                condition_by_sample = networkGroupIDs(object))

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

    message("Estimating model for ", appendLF = FALSE)
    message(k, appendLF = FALSE);message("...", appendLF = TRUE)

    #fit the networks
    fit <- adjDGlasso_minimal(t(data_split_by_condition[[k]]),
                              weights = model_weight_values[[k]],
                              lambda = optimized_lambda)

    #grab the adjacency matrices
    weighted_adjacency_matrices[[k]] <- matrix(data = fit$Theta.glasso,
                                               nrow = num_features, ncol = num_features,
                                               dimnames = list(featureNames(object, original = FALSE),
                                                               featureNames(object, original = FALSE)))
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
#' This function clusters the jointly estimated adjacency matrix constructed using \code{\link{getNetworks}} via
#' the consensus clustering approach described in Ma et al (\emph{Please see the \strong{\emph{Details}} section for more
#' information}) to identify metabolic modules, aka subnetworks, present in the larger networks.
#' Only subnetworks with consensus that meets or exceeds tau are identified as real.
#'
#' @param object A \code{DNEAresults} object
#' @param tau The % agreement threshold among the clustering algorithms for a node to be included in a subnetwork
#' @param max_iterations The maximum number of replicates of the clustering algorithms to perform before consensus is reached
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{createDNEAobject}}, \code{\link{BICtune}}, \code{\link{stabilitySelection}},
#' \code{\link{edgeList}}, \code{\link{adjacencyMatrix}}, \code{\link{getNetworks}}, \code{\link{edgeList}},
#' \code{\link{filterNetworks}}, \code{\link{plotNetworks}}
#'
#' @references
#' Ma J, Karnovsky A, Afshinnia F, Wigginton J, Rader DJ, Natarajan L, Sharma K, Porter AC, Rahman M, He J, Hamm L, Shafi T, Gipson D, Gadegbeku C, Feldman H, Michailidis G, Pennathur S. Differential network enrichment analysis reveals novel lipid pathways in chronic kidney disease. Bioinformatics. 2019 Sep 15;35(18):3441-3452. doi: 10.1093/bioinformatics/btz114. PMID: 30887029; PMCID: PMC6748777. \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6748777/}
#'
#' @details
#' Seven clustering algorithms from the \code{\link{igraph}} package are utilized in this consensus clustering approach:
#' \enumerate{
#' \item \code{\link[igraph:cluster_edge_betweenness]{cluster_edge_betweenness}}
#' \item \code{\link[igraph:cluster_fast_greedy]{cluster_fast_greedy}}
#' \item \code{\link[igraph:cluster_infomap]{cluster_infomap}}
#' \item \code{\link[igraph:cluster_label_prop]{cluster_label_prop}}
#' \item \code{\link[igraph:cluster_louvain]{cluster_louvain}}
#' \item \code{\link[igraph:cluster_walktrap]{cluster_walktrap}}
#' \item \code{\link[igraph:cluster_leading_eigen]{cluster_leading_eigen}}}
#'
#' For each iteration, node membership in each respective cluster is compared across the algorithms, and only the clusters with % agreement
#' greater than tau are kept. A new adjacency graph is then created and clustering is performed again. This occurs iteratively until consensus
#' on stable subnetworks or the specified max_iterations is reached \emph{(Please see references for more details)}.
#'
#' @returns A \code{DNEAresults} object containing sub-network determinations for the nodes within the input network.
#'        A summary of the consensus clustering results can be viewed using \code{\link{CCsummary}}.
#'        Sub-network membership for each node can be found in the "membership" column of the node list, which can be
#'        viewed using \code{\link{nodeList}}.
#'
#' @examples
#' #import completed example data
#' data(TEDDYresults)
#'
#' #identify metabolic modules via consensus clustering
#' TEDDYresults <- clusterNet(object = TEDDYresults, tau = 0.5, max_iterations = 5)
#'
#' #we can also plot the subnetworks
#' plotNetworks(object = TEDDYresults, type = "subnetworks", subtype = 1)
#'
#' @import igraph
#' @export
clusterNet <- function(object,
                       tau = 0.5,
                       max_iterations = 5){

  #test for proper inputs
  if(!inherits(object, "DNEAresults")) stop('the input object should be of class "DNEAresults"!')

  if(tau < 0.5 | tau > 1.0) stop("tau corresponds to a percent agreement among the clustering methods. ",
                                 "As such, tau must be greater than 0.5 and less than 1!",
                                 "Clustering results below this threshold are not reliable -",
                                 "Please see user documentation for more information!")

  if(max_iterations < 1) stop("max_iterations should be a positive integer!")

  #####################################
  #**Join the two condition networks**#
  #####################################

  #create list to hold graph from adjacency matrix
  adjacency_matrix_graphs <- vector("list", length(adjacencyMatrix(object, weighted = TRUE)))
  names(adjacency_matrix_graphs) <- names(adjacencyMatrix(object, weighted = TRUE))

  for (loop_el in names(adjacencyMatrix(object, weighted = TRUE))) {

    adjacency_graph <- graph_from_adjacency_matrix(adjacencyMatrix(object, weighted = TRUE)[[loop_el]], mode = "undirected", weighted = TRUE)
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
  fit <- run_consensus_cluster(joint_graph, tau = tau, max_iterations = max_iterations)
  consensus_membership <- fit$final_consensus_cluster

  #initiate output matrix
  subnetwork_results <- matrix(0, nrow=length(unique(consensus_membership)), numFeatures(object),
                               dimnames = list(paste0("subnetwork", seq(1, length(unique(consensus_membership)))),
                                               vapply(seq(1, length(joint_graph)), function(x) names(joint_graph[[x]]), character(1))))
#####
  #gather results
  for (j in seq(1, nrow(subnetwork_results))){

    #grab features in this subnetwork
    subnetwork_nodes <- consensus_membership == j


    if(sum(subnetwork_nodes) == 1){

      #if subnetwork is only one feature, relabel to "independent"
      consensus_membership[consensus_membership == j] <- "independent"
    }else{

      #concatenate the results
      subnetwork_results[j, consensus_membership == j] <- 1
    }
  }

  #remove empty rows
  subnetwork_results <- subnetwork_results[rowSums(subnetwork_results) != 0, ]

  #update consensus subnetworks
  consensus_membership <- match(consensus_membership, as.numeric(gsub('subnetwork','',rownames(subnetwork_results))))
  consensus_membership[is.na(consensus_membership)] <- "independent"

  #update subnetwork_results rownames
  rownames(subnetwork_results) <- paste0("subnetwork", seq(1, nrow(subnetwork_results)))

  #concatenate the netGSA results table
  summary_list <- list()
  for (loop_cluster in seq(1, nrow(subnetwork_results))){

    cluster_c <- induced.subgraph(joint_graph, V(joint_graph)$name[(subnetwork_results[loop_cluster,] == 1)])
    summary_list[[loop_cluster]] <- data.frame("number_of_nodes"=length(V(cluster_c)),
                                               "number_of_edges"=length(E(cluster_c)),
                                               "number_of_DE.nodes"=sum(as.numeric(table(V(cluster_c)$DE)[names(table(V(cluster_c)$DE))==TRUE])),
                                               "number_of_DE.edges"=sum(as.numeric(table(E(cluster_c)$color)[names(table(E(cluster_c)$color)) %in% c("red", "green")])),
                                               check.names = FALSE)
  }

  summary_stat <- data.frame("subnetworks"= rownames(subnetwork_results), do.call(rbind, summary_list), check.names = FALSE)

  #add independent features
  summary_stat <- rbind(summary_stat,
                        list("independent",
                             sum(consensus_membership == "independent"),
                             0,
                             sum(nodeList(object)$DEstatus[consensus_membership == "independent"]),
                             0))

  #add results to DNEAresults object
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
#' This function performs pathway enrichment analysis on the metabolic modules identified via \code{\link{clusterNet}}
#' using the \code{\link[netgsa:NetGSA]{netgsa::NetGSA()}} algorithm.
#'
#' @param object A \code{DNEAresults}
#' @param min_size The minimum size of metabolic modules for enrichment analysis
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{netGSAresults}}, \code{\link{plotNetworks}}, \code{\link{nodeList}}
#'
#' @references
#' Hellstern M, Ma J, Yue K, Shojaie A. netgsa: Fast computation and interactive visualization for topology-based pathway enrichment analysis. PLoS Comput Biol. 2021 Jun 11;17(6):e1008979. doi: 10.1371/journal.pcbi.1008979. PMID: 34115744; PMCID: PMC8221786. url{https://pubmed.ncbi.nlm.nih.gov/34115744/}
#'
#'
#' @returns A \code{DNEAresults} object after populating the @@netGSA slot. Pathway expression differences for
#'          each node can be found in the node_list. A summary of the NetGSA results can be viewed
#'          using \code{\link{netGSAresults}}.
#'
#' @examples
#' #import completed example data
#' data(TEDDYresults)
#'
#' #perform pathway enrichment analysis using netGSA
#' TEDDYresults <- runNetGSA(object = TEDDYresults, min_size = 5)
#'
#' #view the results
#' netGSAresults(TEDDYresults)
#'
#'
#' @importFrom stats p.adjust
#' @import igraph
#' @importFrom netgsa NetGSA
#' @export
runNetGSA <- function(object, min_size = 5){

  #################################
  #**Prepare data and run netGSA**#
  #################################

  #test for proper input
  if(!inherits(object, "DNEAresults")) stop('the input object should be of class "DNEAresults"!')
  if(min_size <1) stop("min_size parameter should be a positive integer greater than zero!")

  ##set input variables
  adjacency_matrices <- list(list(adjacencyMatrix(x = object, weighted = TRUE)[[1]]),
                             list(adjacencyMatrix(x = object, weighted = TRUE)[[2]]))

  expression_data <- expressionData(object, normalized = FALSE)
  data_groups <- ifelse(networkGroupIDs(object) == networkGroups(object)[1], 1, 2)
  subnetworks <- as.matrix(subnetworkMembership(object))

  ##filter subnetworks to only include those greater than or equal to min_size
  sub_membership <- subnetworkMembership(object)
  filtered_subnetworks <- as.matrix(sub_membership[rowSums(sub_membership) >= min_size, ])

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


  #rename subnetworks
  cluster_names <- CCsummary(object)[-match("independent", CCsummary(object)$subnetworks), ]
  new_cluster_order <- data.frame(subnetworks = c(res$subnetworks, cluster_names$subnetworks[!(cluster_names$subnetworks %in% res$subnetworks)]))

  #update consensus subnetworks
  nodeList(object)$membership <- match(nodeList(object)$membership, as.numeric(gsub('subnetwork','',new_cluster_order$subnetworks)))
  nodeList(object)$membership[is.na(nodeList(object)$membership)] <- "independent"

  #update res subnetwork names
  res$subnetworks <- paste0("subnetwork", seq(1, nrow(res)))

  #update cluster summary subnetwork names
  cluster_names <- CCsummary(object)
  rownames(cluster_names) <- cluster_names$subnetworks
  cluster_names <- cluster_names[new_cluster_order$subnetworks, ]
  rownames(cluster_names) <- NULL
  cluster_names$subnetworks <- paste0("subnetwork", seq(1, nrow(cluster_names)))
  cluster_names <- rbind(cluster_names, CCsummary(object)[match("independent", CCsummary(object)$subnetworks), ])

  #update subnetwork_membership matrix
  sub_membership <- sub_membership[new_cluster_order$subnetworks, ]
  rownames(sub_membership) <- paste0("subnetwork", seq(1, nrow(sub_membership)))

  #update DNEAobject
  netGSAresults(object) <- res
  CCsummary(object) <- cluster_names
  subnetworkMembership(object) <- sub_membership

  #check valid object
  validObject(object)

  return(object)
}
