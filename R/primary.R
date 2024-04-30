#' @include all-methods.R
#' @include all-generics.R
#' @include JSEM-internals.R
#' @include utilities-exported.R
#' @include utilities-internals.R
#' @include clustering-internals.R
#' @include all-classes.R
#' @include start-here.R
NULL

BICtune.DNEAobj <- function(object,
                            lambda_values,
                            interval=1e-3,
                            informed=TRUE,
                            eps_threshold=1e-06,
                            eta_value=0.1,
                            BPPARAM=bpparam(),
                            BPOPTIONS=bpoptions()){


  ##test for proper input
  if(!inherits(object, "DNEAobj")) stop('the input object should be of class "DNEAobj"!')
  if(!is.logical(informed)) stop('"informed" parameter should be TRUE or FALSE!')
  if(interval < 0 | interval > 0.01) stop('"interval" should be between 0 and 0.1!')

  ##initialize input parameters
  dat <- expressionData(x=object, assay="scaled_expression_data")
  n4cov <- max(vapply(dat, ncol, numeric(1)))
  trainX <- t(do.call(cbind, dat))
  trainY <- c(rep(1, ncol(dat[[1]])),
              rep(2, ncol(dat[[2]])))
  bic1 <- NULL
  bic2 <- NULL

  message("Optimizing the lambda hyperparameter using Bayesian-Information Criterion outlined in Guo et al. (2011)\n",
          "A Link to this reference can be found in the function documentation by running ?BICtune() in the console\n")

  ##Pre-define a range of lambda values to evaluate during optimization if none are provided
  if(missing(lambda_values)){

    bic1 <- lambda_tune_dispatch(informed=informed,
                                 FUN='CGM_AHP_tune',
                                 trainX=trainX,
                                 testX=trainX,
                                 trainY=trainY,
                                 BIC=TRUE,
                                 eps=eps_threshold,
                                 eta=eta_value,
                                 asymptotic_lambda=sqrt(log(numFeatures(object))/n4cov),
                                 interval=interval,
                                 BPPARAM=BPPARAM,
                                 BPOPTIONS=BPOPTIONS)

    lambda_values <- bic1$new_lambda_values
    bic1 <- bic1$bic
  }else{

    message("Provided lambda values will be used for optimization...", appendLF=TRUE)
    lambda_values <- unlist(lambda_values)
    if(any(lambda_values < 0 | lambda_values > 1)){

      stop("The lambda parameter should be a value between 0 and 1 only!")
    }
  }

  ##search the estimated area for the optimized lambda
  message("Fine-tuning Lambda...")
  bic2 <- tune_lambda(lambda_values=lambda_values,
                      FUN='CGM_AHP_tune',
                      trainX=trainX,
                      testX=trainX,
                      trainY=trainY,
                      BIC=TRUE,
                      eps=eps_threshold,
                      eta=eta_value,
                      BPPARAM=BPPARAM,
                      BPOPTIONS=BPOPTIONS)

  ##concatenate tuning runs and reorder by lambda value
  lambda_tested <- append(bic1[["lambda_values"]], bic2[["lambda_values"]])
  BIC_guo <- append(bic1[["BIC_guo"]], bic2[["BIC_guo"]])

  ##update DNEA object
  BICscores(object) <- BIC_guo
  optimizedLambda(object) <- bic2[["lastar_guo"]]
  lambdas2Test(object) <- lambda_tested

  message("The optimal Lambda hyper-parameter has been set to: ",
          bic2[["lastar_guo"]], "!")

  #check valid object
  validObject(object)
  return(object)
}

BICtune.matrix <- function(object,
                           lambda_values,
                           interval=1e-3,
                           informed=TRUE,
                           eps_threshold=1e-06,
                           eta_value=0.1,
                           BPPARAM=bpparam(),
                           BPOPTIONS=bpoptions()){


  ##test for proper input
  if(!inherits(object, "matrix")) stop('object should be a list of matrices!')
  if(!is.logical(informed)) stop('"informed" parameter should be TRUE or FALSE!')
  if(interval < 0 | interval > 0.01) stop('"interval" should be between 0 and 0.1!')

  ##initialize input parameters
  n4cov <- ncol(object)
  trainY <- rep(1, n4cov)
  num_features <- ncol(object)
  object <- t(object)
  bic1 <- NULL
  bic2 <- NULL

  ##Pre-define a range of lambda values to evaluate during optimization if none are provided
  if(missing(lambda_values)){

    bic1 <- lambda_tune_dispatch(informed=informed,
                                 FUN='CGM_AHP_tune',
                                 trainX=object,
                                 testX=object,
                                 trainY=trainY,
                                 BIC=TRUE,
                                 eps=eps_threshold,
                                 eta=eta_value,
                                 asymptotic_lambda=sqrt(log(num_features)/n4cov),
                                 interval=interval,
                                 BPPARAM=BPPARAM,
                                 BPOPTIONS=BPOPTIONS)

    lambda_values <- bic1$new_lambda_values
    bic1 <- bic1$bic
  }else{

    message("Provided lambda values will be used for optimization...", appendLF=TRUE)
    lambda_values <- unlist(lambda_values)
    if(any(lambda_values < 0 | lambda_values > 1)){

      stop("The lambda parameter should be a value between 0 and 1 only!")
    }
  }

  ##search the estimated area for the optimized lambda
  message("Fine-tuning Lambda...")
  bic2 <- tune_lambda(lambda_values=lambda_values,
                      FUN='CGM_AHP_tune',
                      trainX=object,
                      testX=object,
                      trainY=trainY,
                      BIC=TRUE,
                      eps=eps_threshold,
                      eta=eta_value,
                      BPPARAM=BPPARAM,
                      BPOPTIONS=BPOPTIONS)

  ##concatenate tuning runs and reorder by lambda value
  lambda_tested <- append(bic1[["lambda_values"]], bic2[["lambda_values"]])
  BIC_guo <- append(bic1[["BIC_guo"]], bic2[["BIC_guo"]])

  message("The optimal Lambda hyper-parameter has been set to: ",
          bic2[["lastar_guo"]], "!")

  return(list(BIC_guo=BIC_guo,
              optimized_lambda=bic2[["lastar_guo"]],
              tested_lambda_values=lambda_tested))
}

#' Optimize the lambda regularization parameter for the glasso-based
#' network models using Bayesian-information Criterion
#'
#' This function will calculate the Bayesian information criterion (BIC)
#' and likelihood for a range of lambda values that are automatically
#' generated (\emph{please see \strong{Details} for more info}) or that are
#' user-specified. The lambda value with the minimum BIC score is the optimal
#' lambda value for the data set and is stored in the DNEAobj object for use in
#' stability selection using \code{\link{stabilitySelection}}.
#'
#' @param object A \code{\link{DNEAobj}} object. See \code{\link{createDNEAobject}}
#' @param lambda_values \emph{\strong{OPTIONAL -}} A list of values to test while optimizing
#' the lambda parameter. If not provided, a set of lambda values are chosen
#' based on the theoretical value for the asymptotically valid lambda. More
#' information about this can be found in the details section.
#'
#' @param interval A numeric value indicating the specificity by which to
#' optimize lambda. The default value is 1e-3, which indicates lambda will
#' be optimized to 3 decimal places. The value should be between 0 and 0.1.
#'
#' @param informed TRUE/FALSE indicating whether the asymptotic properties
#' of lambda for large data sets should be utilized to tune the parameter.
#' This reduces the necessary number of computations for optimization.

#' @param eps_threshold A significance cut-off for thresholding network
#' edges. The default value is 1e-06.
#' This value generally should not change.
#'
#' @param eta_value A tuning parameter that that ensures that the empirical
#' covariance matrix of the data is positive definite so that we can
#' calculate its inverse. The default value is 0.01.
#'
#' @param BPPARAM A \code{\link{BiocParallel}} object.
#'
#' @param BPOPTIONS a list of options for BiocParallel created using
#' the \code{\link[BiocParallel:bpoptions]{bpoptions}} function.
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{optimizedLambda}},
#' \code{\link[BiocParallel:bpparam]{bpparam}},
#' \code{\link[BiocParallel:bpoptions]{bpoptions}}
#' \code{\link[glasso:glasso]{glasso}}
#'
#' @references Guo J, Levina E, Michailidis G, Zhu J.
#' Joint estimation of multiple graphical models.
#' Biometrika. 2011 Mar;98(1):1-15. doi: 10.1093/biomet/asq060.
#' Epub 2011 Feb 9. PMID: 23049124; PMCID: PMC3412604.
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3412604/}
#'
#' @details
#' There are several ways to optimize the lambda parameter for a
#' glasso model - We utilize Bayesian-information criterion (BIC) to
#' optimize the lambda parameter in DNEA because it is a more balanced
#' method and less computationally expensive. We can reduce the total
#' number of values that need to be tested in optimization by carefully
#' selecting values around the asymptotically valid lambda for data
#' sets with many samples and many features following the equation:
#'  \deqn{\lambda = \sqrt{ \ln (num. features) / num. samples}}{ lambda = sqrt(ln(num. features) / num. samples)}
#'
#' For smaller data sets, the asymptotically valid lambda is described
#' by modifying the previous equation to include an unknown constant, c,
#' that needs to be determined mathematically. Therefore, to optimize
#' lambda we modify the previous equation as follows:
#' \deqn{\lambda = c \sqrt{ \ln (num. features) / num. samples}}{lambda = c*sqrt(ln(num. features) / num. samples)}
#'
#' where c takes on values between 0 and the
#' theoretical maximum of C in intervals of 0.02. C is then estimated
#' and a new range is tested to the specificity of the "interval"
#' input. More information regarding the optimization method
#' deployed here can be found in the Guo et al. (2011)
#' paper referenced below.
#'
#' @returns A \code{\link{DNEAobj}} object containing the BIC and likelihood
#' scores for every lambda value tested, as well as the
#' optimized lambda value
#'
#' @examples
#' #import BiocParallel package
#' library(BiocParallel)
#'
#' #load example data
#' data(TEDDY)
#' data(T1Dmeta)
#'
#' #make sure metadata and expression data are in same order
#' T1Dmeta <- T1Dmeta[colnames(TEDDY),]
#'
#' #create group labels
#' group_labels <- T1Dmeta$group
#' names(group_labels) <- rownames(T1Dmeta)
#'
#' #initiate DNEAobj
#' dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
#'                             group_labels = group_labels)
#'
#' #optimize lambda parameter
#' dnw <- BICtune(object=dnw,
#'                informed=TRUE,
#'                interval=0.01)
#'
#' @import glasso
#' @importFrom BiocParallel bplapply bpparam bpoptions bptasks
#' @rdname BICtune-methods
#' @aliases BICtune
#' @export
setMethod("BICtune", signature(object="DNEAobj"), BICtune.DNEAobj)

#' @rdname BICtune-methods
#' @aliases BICtune
#' @export
setMethod("BICtune", signature(object="matrix"), BICtune.matrix)

#' Stability selection calculates selection probabilities for every
#' possible feature-feature interaction within the input data
#'
#' This function randomly samples the input data and fits a glasso model
#' with the sampled data for \strong{\emph{nreps}} number of replicates. The
#' resulting adjacency matrices are summed together and selection probabilities
#' for each feature-feature interaction are calculated. Stability selection is
#' particularly useful for smaller data sets. A large number of replicates
#' should be performed (the default is 1000). The exact method deployed
#' varies slightly whether or not additional sub-sampling of the data is
#' performed. More information can be found in the
#' \strong{\emph{Details}} section.
#'
#' @param object A \code{\link{DNEAobj}} object.
#'
#' @param subSample TRUE/FALSE indicating whether the number of samples
#' are unevenly split by condition and subsampling should be performed
#' when randomly sampling to even out the groups.
#'
#' @param nreps The total number of replicates to perform in stability
#' selection. The default is 1000.
#'
#' @param optimal_lambda \emph{\strong{OPTIONAL -}} The optimal lambda value to be
#' used in the model. This parameter is only necessary if
#' \code{\link{BICtune}} is not performed prior.
#'
#' @param BPPARAM a BiocParallel object.
#'
#' @param BPOPTIONS a list of options for BiocParallel created using
#' the \code{\link[BiocParallel:bpoptions]{bpoptions}} function.
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{selectionProbabilities}},
#' \code{\link{selectionResults}},
#' \code{\link[BiocParallel:bpparam]{bpparam}},
#' \code{\link[BiocParallel:bpoptions]{bpoptions}}
#' \code{\link[glasso:glasso]{glasso}}
#'
#'
#' @references Ma J, Karnovsky A, Afshinnia F, Wigginton J, Rader DJ,
#' Natarajan L, Sharma K, Porter AC, Rahman M, He J, Hamm L, Shafi T,
#' Gipson D, Gadegbeku C, Feldman H, Michailidis G, Pennathur S.
#' Differential network enrichment analysis reveals novel lipid
#' pathways in chronic kidney disease. Bioinformatics.
#' 2019 Sep 15;35(18):3441-3452. doi: 10.1093/bioinformatics/btz114.
#' PMID: 30887029; PMCID: PMC6748777.
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6748777/}
#'
#' Nicolai, M., & Peter, B. (2010). Stability selection.
#' Journal of the Royal Statistical Society: Series B
#' (Statistical Methodology), 72(4), 417-473.
#' \url{https://stat.ethz.ch/Manuscripts/buhlmann/stability.pdf}
#'
#' @details
#' Stability selection provides an additional approach by which to
#' regularize the network model and create more robust results, particularly
#' when \strong{\emph{p >> n}}. Stability selection works by randomly
#' sampling (without replacement) the input data many times and fitting
#' a glasso model to each subset of sampled data. The unwieghted adjacency
#' matrix from each model is summed together (A feature-feature interaction
#' is considered present if the partial correlation value is above 1e-5),
#' and the probability of an edge being selected in a random subset of the
#' data is calculated by dividing the number of times an edge was selected in
#' the replicates over the total number of replicates. This results in a
#' selection probability for every possible feature-feature interaction
#' that is used to modify the regularization parameter
#' via the following equation:
#' \deqn{\rho = \lambda*(1 / (0.0001 + selection.probability))}{ rho = lambda*(1 / (0.0001 + selection.probability))}
#'
#'
#' However, when the sample groups are very unbalanced, randomly
#' sampling strongly favors the larger group, resulting in over
#' representation. In order to combat this, setting subSample=TRUE
#' modifies the random sample by sub-sampling the experimental groups
#' to even out the sample numbers. In this method, 90% of the smaller
#' group is randomly sampled without replacement, and an
#' additional 10% is randomly sampled without replacement from
#' the entire group to preserve the variance. The larger group
#' is randomly sampled to have 1.3 times the number of samples
#' present in the smaller group. This method ensures that each
#' group is equally represented in stability selection. \cr
#'
#' The principles of stability selection remain similar with both methods,
#' however, there are a few small differences. Stability selection
#' \emph{without} additional sub-sampling randomly samples 50% of each group
#' (without replacement) and fits a model for both halves of the sampled data.
#' Since nearly all of the data for the smaller group is used \emph{with}
#' additional sub-sampling, only one model is fit per replicate when
#' subSample=TRUE. This means that at the default value of nreps=500,
#' 1000 randomly sampled models are fit in total \emph{without} sub-sampling,
#' but 500 randomly sampled models are fit in total \emph{with} sub-sampling.
#' More details about the stability approach deployed in this function can be
#' found in Ma et al. (2019) referenced below.
#'
#'
#' @returns A \code{\link{DNEAobj}} object after populating the
#' stable_networks slot of the object. It contains the selection
#' results from stability selection as well as the calculated
#' selection probabilities.
#'
#' @examples
#' #import BiocParallel package
#' library(BiocParallel)
#'
#' #dnw is a DNEAobj with the results generated for the example data
#' #accessed by running data(TEDDY) in the console. The workflow
#' #for this data can be found in the vignette accessed by
#' #running browseVignettes("DNEA") in the console.
#' data(dnw)
#'
#' # perform stability selection
#' dnw <- stabilitySelection(object=dnw,
#'                           subSample=FALSE,
#'                           nreps=4,
#'                           BPPARAM=bpparam())
#'
#' @import glasso
#' @importFrom BiocParallel bplapply bpparam bpoptions bptasks
#' @export
stabilitySelection <- function(object,
                               subSample=FALSE,
                               nreps=500,
                               optimal_lambda,
                               BPPARAM=bpparam(),
                               BPOPTIONS=bpoptions()){

  ##test for proper input
  if(!inherits(object, "DNEAobj")){
    stop('the input object should be of class "DNEAobj"!')
  }
  if(nreps < 1 | !is.numeric(nreps)) {
    stop("nreps must be positive - ",
    "We recommend setting nreps to 500 at minimum!")
  }
  if(!is.logical(subSample)) {
    stop('"subSample" parameter should be TRUE or FALSE!')
  }

  ## stabilitySelection requires lambda hyper-parameter. Will use
  ## optimal_lambda if supplied, otherwise looks for
  ## @hyperparameter[["optimized_lambda"]] in DNEAobject
  if(!missing(optimal_lambda)){

    if(optimal_lambda < 0 | optimal_lambda > 1){

      stop("The lambda parameter should be a value between 0 and 1 only!")
    }
    ##set lambda for use downstream
    optimized_lambda <- optimal_lambda
    message('IMPORTANT: optimal_lambda argument was provided - The default parameters will ',
            'be overrided and optimal_lambda will be used in analysis')

    if(is.null(optimizedLambda(object))){
      optimizedLambda(object) <- optimal_lambda
      message('@hyperparameter[["optimized_lambda"]] was previously empty and now set to optimal_lambda argument')
    }
  }else if(!is.null(optimizedLambda(object))){
    optimized_lambda <- optimizedLambda(object)
    message('The lambda value stored in the DNEAobj will be used for analysis (this can be ',
            'accessed via the optimizedLambda() function')
  }else{
    ##setting optimized_lambda=NULL will default to a lambda of
    ##sqrt(log(# features) / # samples) in adjDGlasso_minimal
    optimized_lambda <- NULL

    warning("No lambda value was supplied for the model - sqrt(log(# features) / # samples) will be used in the analyis. ",
            "However, We highly recommend optimizing the lambda parameter by running BICtune(), ",
            "or providing a calibrated lambda value using the optimal_lambda parameter prior to analysis ",
            "if the dataset contains ~500 or more samples.")
  }

  data_split_by_condition <- lapply(expressionData(x=object, assay="scaled_expression_data"),
                                    function(d) t(d))

  message("Using Lambda hyper-parameter: ", optimized_lambda, "!\n",
          "stabilitySelection will be performed with ", nreps, " replicates!")

  ##set internal function to use
  if(subSample){

    message("Additional sub-sampling will be performed on uneven groups")
    ss_function <- "CGM_AHP_stabsel_subsample"
  }else if(!subSample){

    message("No additional sub-sampling will be performed. Sample groups will both be randomly sampled 50%")
    ss_function <- "CGM_AHP_stabsel"
  }

  ##initialize parameters and run stability selection
  stabsel_init_param <- stabsel_init(listX=data_split_by_condition, nreps=nreps)
  stab_sel <- BiocParallel:: bplapply(X=seq(1, nreps),
                                      FUN=ss_function,
                                      init_param=stabsel_init_param,
                                      listX=data_split_by_condition,
                                      lastar=optimized_lambda,
                                      BPPARAM=BPPARAM,
                                      BPOPTIONS=BPOPTIONS)

  ##concatenate results for output
  selection_results <- vector("list", length(networkGroups(object)))
  names(selection_results) <- networkGroups(object)

  selection_probabilities <- vector("list", length(networkGroups(object)))
  names(selection_probabilities) <- networkGroups(object)

  for (k in seq(1, length(selection_results))){
    selection_results[[k]] <- lapply(stab_sel, function(r) r$mat[[k]])
    selection_results[[k]] <- Reduce("+", selection_results[[k]])

    if (subSample){

      message("Calculating selection probabilities WITH subsampling for...",
              names(selection_results)[[k]], "...", appendLF=TRUE)
      selection_probabilities[[k]] <- selection_results[[k]]/(nreps)
    } else {

      message("Calculating selection probabilities WITHOUT subsampling for...",
              names(selection_results)[[k]], "...", appendLF=TRUE)
      selection_probabilities[[k]] <- selection_results[[k]]/(2 * nreps)
    }
  }
  selectionResults(object) <- selection_results
  selectionProbabilities(object) <- selection_probabilities

  ##check valid object
  validObject(object)
  return(object)
}

#' Construct the GLASSO-based biological Networks
#'
#' This function constructs a biological network for each experimental
#' condition using the joint estimation method described
#' in Ma et al. (2019) (\emph{please see references below}). If
#' \code{\link{stabilitySelection}} was performed previously,
#' the selection probabilities will be used to for model optimization
#' when constructing the networks (please see the \strong{\emph{Details}}
#' section of \code{\link{stabilitySelection}} for more information).
#'
#'
#' @param object A \code{\link{DNEAobj}} object.
#'
#' @param lambda_values **OPTIONAL** A list of values to test while optimizing
#' the lambda parameter. If not provided, a set of lambda values are chosen
#' based on the theoretical value for the asymptotically valid lambda. More
#' information about this can be found in the details section of
#' \code{\link{BICtune}}.
#'
#' @param aprox TRUE/FALSE indicating whether \code{\link{BICtune}}
#' should be used to optimize the lambda value for each condition. If
#' `aprox==FALSE`, sqrt(log(# features)/#samples) is used to approximate
#' lambda.
#'
#' @param optimal_lambdas \emph{\strong{OPTIONAL -}} The lambda value
#' to be used in analysis. If not provided, the lambda value is
#' determined based on the input of the "aprox" parameter.

#' @param informed TRUE/FALSE indicating whether the asymptotic properties
#' of lambda for large data sets should be utilized to tune the parameter.
#' This reduces the necessary number of computations for optimization.
#'
#' @param interval A numeric value indicating the specifity by which to
#' optimize lambda. The default value is 1e-3, which indicates lambda will
#' be optimized to 3 decimal places. The value should be between 0 and 0.1.
#'
#' @param eps_threshold A significance cut-off for thresholding network edges.
#'        The default value is 1e-06. This value generally should not change.
#'
#' @param eta_value A tuning parameter that that ensures that the empirical
#' covariance matrix of the data is positive definite so that we can
#' calculate its inverse. The default value is 0.01.
#'
#' @param BPPARAM A \code{\link{BiocParallel}} object.
#'
#' @param BPOPTIONS a list of options for BiocParallel created using
#' the \code{\link[BiocParallel:bpoptions]{bpoptions}} function.
#'
#' @param eps_threshold A numeric value between 0 and 1 by which to
#' threshold the partial correlation values for edge identification.
#' Edges with an absolute partial correlation value below this
#' threshold will be zero'd out from the adjacency matrix.
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{edgeList}},\code{\link{adjacencyMatrix}}
#'
#' @references
#' Ma J, Karnovsky A, Afshinnia F, Wigginton J, Rader DJ, Natarajan L,
#' Sharma K, Porter AC, Rahman M, He J, Hamm L, Shafi T, Gipson D,
#' Gadegbeku C, Feldman H, Michailidis G, Pennathur S.
#' Differential network enrichment analysis reveals novel lipid pathways in
#' chronic kidney disease. Bioinformatics. 2019 Sep 15;35(18):3441-3452.
#' doi: 10.1093/bioinformatics/btz114. PMID: 30887029; PMCID: PMC6748777.
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6748777/}
#'
#' Iyer GR, Wigginton J, Duren W, LaBarre JL, Brandenburg M, Burant C,
#' Michailidis G, Karnovsky A. Application of Differential Network
#' Enrichment Analysis for Deciphering Metabolic Alterations. Metabolites.
#' 2020 Nov 24;10(12):479. doi: 10.3390/metabo10120479. PMID: 33255384;
#' PMCID: PMC7761243. \url{https://pubmed.ncbi.nlm.nih.gov/33255384/}
#'
#' @returns A \code{\link{DNEAobj}} object after populating the adjaceny_matrix
#' and edge_list slots with the corresponding adjacency_matrix for each
#' sample condition as well as the network edge list.
#'
#' @examples
#' #dnw is a DNEAobj with the results generated for the example data
#' #accessed by running data(TEDDY) in the console. The workflow
#' #for this data can be found in the vignette accessed by
#' #running browseVignettes("DNEA") in the console.
#' data(dnw)
#'
#' #construct the networks
#' dnw <- getNetworks(object=dnw, aprox = TRUE)
#'
#' #now we can plot the group networks
#' plotNetworks(object=dnw, type="group_networks")
#'
#' @importFrom gdata lowerTriangle
#' @importFrom utils combn
#' @export
#'
getNetworks <- function(object,
                        lambda_values,
                        aprox=FALSE,
                        informed=TRUE,
                        interval=1e-3,
                        eps_threshold=1e-06,
                        eta_value=0.1,
                        optimal_lambdas,
                        BPPARAM=bpparam(),
                        BPOPTIONS=bpoptions()){

  ##test for proper input
  if(!inherits(object, "DNEAobj")){
    stop('the input object should be of class "DNEAobj"!')
  }
  if(eps_threshold <=0 | eps_threshold >= 1) {
    stop("eps_threshold should be between 0 and 1 only!")
  }

  ##initialize input parameters
  data_split_by_condition <- expressionData(x=object, assay="scaled_expression_data")
  num_samples <- vapply(data_split_by_condition, function(x) ncol(x), FUN.VALUE=numeric(1))
  names(num_samples) <- names(data_split_by_condition)
  num_features <- numFeatures(object)
  Ip <- diag(rep(1, num_features))

  ##initiate output data structures
  tuned_lambdas <- vector("list", length(networkGroups(object)))
  names(tuned_lambdas) <- networkGroups(object)

  weighted_adjacency_matrices <- vector("list", length(networkGroups(object)))
  names(weighted_adjacency_matrices) <- networkGroups(object)

  unweighted_adjacency_matrices <- vector("list", length(weighted_adjacency_matrices))
  names(unweighted_adjacency_matrices) <- names(weighted_adjacency_matrices)

  if(!missing(optimal_lambdas)){

    if(aprox) warning("Lambda values were provided even though aprox=TRUE. ",
                      "The provided lambda values will be used in analysis. ",
                      "To aproximate lambda, remove the supplied values.")

    optimized_lambdas <- optimal_lambdas[names(data_split_by_condition)]
    if(length(optimal_lambdas != 2)){
      stop("Two lambda values should be supplied - one for each condition!")
    }
    if(optimal_lambdas < 0 | optimal_lambdas > 1){
      stop("The lambda parameter should be a value between 0 and 1 only!")
    }

  }else if(missing(optimal_lambdas)){
    if(aprox){

      message("Lambda will be aproximated by sqrt(log(# features)/# samples)\n",
              "NOTE: If your dataset contains 500 or more samples per ",
              "experimental condition, you should set ",
              ' "aprox=FALSE" and tune lambda for each network!')

      optimized_lambdas <- sqrt(log(num_features) / num_samples)

    }else{
      message("Optimizing the lambda hyperparameter using Bayesian-Information Criterion outlined in Guo et al. (2011)\n",
              "A Link to this reference can be found in the function documentation by running ?BICtune() in the console\n",
              "NOTE: if your dataset contains fewer than ~500 samples per experimental condition, consider setting ",
              ' "aprox=TRUE". This will provide more reliable results')

      optimized_lambdas <- NULL
      for(x in names(data_split_by_condition)){

        message("\nTUNING LAMBDA FOR ", x,"!:\n",
                "--------------------------------------------------\n")
        tuned_lambdas[[x]] <- BICtune(object=data_split_by_condition[[x]],
                                      informed=informed,
                                      interval=interval,
                                      eps_threshold=eps_threshold,
                                      eta_value=eta_value,
                                      BPPARAM=BPPARAM,
                                      BPOPTIONS=BPOPTIONS)
        optimized_lambdas[[x]] <- tuned_lambdas[[x]][["optimized_lambda"]]
      }
    }
  }

  ##model will use selection weights if provided
  if (!is.null(selectionProbabilities(object))){

    ##grab selection probabilities and calculate rho
    selection_prob <- selectionProbabilities(object)
    model_weight_values <- vector(mode="list", length=2)
    for(x in seq(1, length(model_weight_values))){

      model_weight_values[[x]] <- 1/(1e-04 + as.matrix(selection_prob[[x]]))
    }
    names(model_weight_values) <- names(selection_prob)
    message('selection_probabilites from stability selection will be used in glasso model!\n')
  } else{

    message("No selection_probabilities were found. We recommend running
            stabilitySelection() prior to estimating the glasso model!\n")
    model_weight_values <- list(matrix(rep(1, num_features^2), num_features, num_features),
                                matrix(rep(1, num_features^2), num_features, num_features))
  }

  names(model_weight_values) <- names(selectionProbabilities(object))

  ##estimate the partial correlation matrix for each condition
  for (k in networkGroups(object)){

    message("Estimating model for ", k, "...using ",
            optimized_lambdas[[k]], " for lambda...")

    fit <- adjDGlasso_minimal(t(data_split_by_condition[[k]]),
                              weights=model_weight_values[[k]],
                              lambda=optimized_lambdas[[k]])

    weighted_adjacency_matrices[[k]] <- matrix(data=fit$Theta.glasso,
                                               nrow=num_features, ncol=num_features,
                                               dimnames=list(featureNames(object, original=FALSE),
                                                             featureNames(object, original=FALSE)))
  }

  ##input weighted_adjacency_matrices into object
  adjacencyMatrix(x=object, weighted=TRUE) <- weighted_adjacency_matrices

  ## filter the weighted_adjacency_matrices by eps_threshold
  ## and create unweighted_adjacency_matrices plus edge list
  object <- filterNetworks(data=object, pcor=eps_threshold)

  #check valid object
  validObject(object)
  return(object)
}

#' Identify metabolic modules within the biological networks using a
#' consensus clustering approach
#'
#' This function clusters the jointly estimated adjacency matrices
#' constructed using \code{\link{getNetworks}} via the consensus clustering
#' approach described in Ma et al. (\emph{Please see the
#' \strong{\emph{Details}} section for more information}) to identify
#' metabolic modules, aka sub networks, present in the larger networks.
#' Only sub networks with consensus that meets or exceeds tau are
#' identified as real.
#'
#' @param object A \code{\link{DNEAobj}} object.
#'
#' @param tau The % agreement among the clustering algorithms
#' for a node to be included in a sub network.
#'
#' @param max_iterations The maximum number of replicates of
#' consensus clustering to be performed if consensus is not
#' reached.
#'
#' @param verbose TRUE/FALSE whether a progress bar should be
#' displayed in the console.
#'
#' @author Christopher Patsalis
#'
#' @seealso \code{\link{plotNetworks}}
#'
#' @references
#' Ma J, Karnovsky A, Afshinnia F, Wigginton J, Rader DJ, Natarajan L,
#' Sharma K, Porter AC, Rahman M, He J, Hamm L, Shafi T, Gipson D,
#' Gadegbeku C, Feldman H, Michailidis G, Pennathur S.
#' Differential network enrichment analysis reveals novel lipid pathways
#' in chronic kidney disease. Bioinformatics. 2019 Sep 15;35(18):3441-3452.
#' doi: 10.1093/bioinformatics/btz114. PMID: 30887029; PMCID: PMC6748777.
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6748777/}
#'
#' @details
#' Seven clustering algorithms from the \code{\link{igraph}} package
#' are utilized in this consensus clustering approach:
#' \enumerate{
#' \item \code{\link[igraph:cluster_edge_betweenness]{cluster_edge_betweenness}}
#' \item \code{\link[igraph:cluster_fast_greedy]{cluster_fast_greedy}}
#' \item \code{\link[igraph:cluster_infomap]{cluster_infomap}}
#' \item \code{\link[igraph:cluster_label_prop]{cluster_label_prop}}
#' \item \code{\link[igraph:cluster_louvain]{cluster_louvain}}
#' \item \code{\link[igraph:cluster_walktrap]{cluster_walktrap}}
#' \item \code{\link[igraph:cluster_leading_eigen]{cluster_leading_eigen}}}
#'
#' For each iteration, node membership in a respective cluster is
#' compared across the algorithms, and only the nodes with tau %
#' agreement for a given cluster are kept. A new adjacency graph is
#' then created and clustering is performed again. This occurs iteratively
#' until consensus on is reached stable sub networks or the specified
#' "max_iterations" is reached
#' \emph{(Please see references for more details)}.
#'
#' @returns A \code{\link{DNEAobj}} object containing sub network
#' determinations for the nodes within the input network. A summary of the
#' consensus clustering results can be viewed using \code{\link{CCsummary}}.
#' Sub network membership for each node can be found in the "membership"
#' column of the node list, which can be accessed using \code{\link{nodeList}}.
#'
#' @examples
#' #dnw is a DNEAobj with the results generated for the example data
#' #accessed by running data(TEDDY) in the console. The workflow
#' #for this data can be found in the vignette accessed by
#' #running browseVignettes("DNEA") in the console.
#' data(dnw)
#'
#' #identify metabolic modules via consensus clustering
#' dnw <- clusterNet(object=dnw, tau=0.5, max_iterations=5)
#'
#' #we can also plot the subnetworks
#' plotNetworks(object=dnw, type="sub_networks", subtype=1)
#'
#' @import igraph
#' @export
clusterNet <- function(object,
                       tau=0.5,
                       max_iterations=5,
                       verbose=TRUE){

  ##test for proper inputs
  if(!inherits(object, "DNEAobj")) {
    stop('the input object should be of class "DNEAobj"!')
  }
  if(tau < 0.5 | tau > 1.0) {
    stop("tau corresponds to a percent agreement among the clustering methods. ",
                                 "As such, tau must be greater than 0.5 and less than 1!",
                                 "Clustering results below this threshold are not reliable -",
                                 "Please see user documentation for more information!")
  }
  if(max_iterations < 1) {
    stop("max_iterations should be a positive integer!")
  }
  if(!is.logical(verbose)) {
    stop('"verbose" parameter should be TRUE or FALSE!')
  }

  ##create list to hold graph from adjacency matrix
  adjacency_matrix_graphs <- vector("list", length(adjacencyMatrix(object, weighted=TRUE)))
  names(adjacency_matrix_graphs) <- names(adjacencyMatrix(object, weighted=TRUE))

  for (loop_el in names(adjacencyMatrix(object, weighted=TRUE))) {

    adjacency_graph <- graph_from_adjacency_matrix(adjacencyMatrix(object, weighted=TRUE)[[loop_el]],
                                                   mode="undirected", weighted=TRUE)
    V(adjacency_graph)$name <- as.character(featureNames(object))
    adjacency_matrix_graphs[[loop_el]] <- adjacency_graph
  }

  ##join adjacency matrix graphs and modify
  joint_graph <- igraph::union(adjacency_matrix_graphs[[1]], adjacency_matrix_graphs[[2]])
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

  ##run consensus cluster
  fit <- run_consensus_cluster(joint_graph, tau=tau, max_iterations=max_iterations, verbose=verbose)
  consensus_membership <- fit$final_consensus_cluster

  ##create output
  subnetwork_results <- matrix(0, nrow=length(unique(consensus_membership)), numFeatures(object),
                               dimnames=list(paste0("subnetwork", seq(1, length(unique(consensus_membership)))),
                                             vapply(seq(1, length(joint_graph)), function(x) names(joint_graph[[x]]), character(1))))
  for (j in seq(1, nrow(subnetwork_results))){

    subnetwork_nodes <- consensus_membership == j
    if(sum(subnetwork_nodes) == 1){

      ##if subnetwork is only one feature, relabel to "independent"
      consensus_membership[consensus_membership == j] <- "independent"
    }else{

      ##concatenate the results
      subnetwork_results[j, consensus_membership == j] <- 1
    }
  }

  ##remove empty rows
  subnetwork_results <- subnetwork_results[rowSums(subnetwork_results) != 0, ]

  ##update consensus subnetworks
  consensus_membership <- match(consensus_membership,
                                as.numeric(gsub('subnetwork','',rownames(subnetwork_results))))
  consensus_membership[is.na(consensus_membership)] <- "independent"

  #update subnetwork_results rownames
  rownames(subnetwork_results) <- paste0("subnetwork", seq(1, nrow(subnetwork_results)))

  #create output table
  summary_list <- list()
  for (loop_cluster in seq(1, nrow(subnetwork_results))){

    cluster_c <- induced.subgraph(joint_graph, V(joint_graph)$name[(subnetwork_results[loop_cluster,] == 1)])
    summary_list[[loop_cluster]] <- data.frame("number_of_nodes"=length(V(cluster_c)),
                                               "number_of_edges"=length(E(cluster_c)),
                                               "number_of_DE.nodes"=sum(as.numeric(table(V(cluster_c)$DE)[names(table(V(cluster_c)$DE))==TRUE])),
                                               "number_of_DE.edges"=sum(as.numeric(table(E(cluster_c)$color)[names(table(E(cluster_c)$color)) %in% c("red", "green")])),
                                               check.names=FALSE)
  }
  summary_stat <- data.frame("subnetworks"= rownames(subnetwork_results),
                             do.call(rbind, summary_list),
                             check.names=FALSE)
  summary_stat <- rbind(summary_stat,
                        list("independent",
                             sum(consensus_membership == "independent"),
                             0,
                             sum(nodeList(object)$DEstatus[consensus_membership == "independent"]),
                             0))

  ##add results to DNEAobj object
  nodeList(object)[["membership"]] <- consensus_membership
  object@consensus_clustering <- new(Class="consensusClusteringResults",
                                     summary=summary_stat,
                                     subnetwork_membership=data.frame(subnetwork_results),
                                     adjacency_graphs=append(adjacency_matrix_graphs,
                                                             list(joint_graph=joint_graph)))
  ##check valid object
  validObject(object)
  return(object)
}

#' Identify metabolic modules that are enriched across
#' experimental conditions using NetGSA
#'
#' This function performs pathway enrichment analysis on the metabolic
#' modules identified via \code{\link{clusterNet}} using the
#' \code{\link[netgsa:NetGSA]{netgsa::NetGSA()}} algorithm.
#'
#' @param object A \code{\link{DNEAobj}}.
#'
#' @param min_size The minimum size of a given metabolic
#' module for to be tested for enrichment across the
#' experimental condition.
#'
#' @author Christopher Patsalis
#'
#' @seealso
#' \code{\link{netGSAresults}}
#' \code{\link{clusterNet}}
#' \code{\link[netgsa:NetGSA]{NetGSA}}
#'
#' @references
#' Hellstern M, Ma J, Yue K, Shojaie A.
#' netgsa: Fast computation and interactive visualization for
#' topology-based pathway enrichment analysis.
#' PLoS Comput Biol. 2021 Jun 11;17(6):e1008979.
#' doi: 10.1371/journal.pcbi.1008979. PMID: 34115744;
#' PMCID: PMC8221786. url{https://pubmed.ncbi.nlm.nih.gov/34115744/}
#'
#'
#' @returns A \code{\link{DNEAobj}} object after populating the @@netGSA
#' slot. A summary of the NetGSA results can be viewed
#' using \code{\link{netGSAresults}}.
#'
#' @examples
#' #dnw is a DNEAobj with the results generated for the example data
#' #accessed by running data(TEDDY) in the console. The workflow
#' #for this data can be found in the vignette accessed by
#' #running browseVignettes("DNEA") in the console.
#' data(dnw)
#'
#' #perform pathway enrichment analysis using netGSA
#' dnw <- runNetGSA(object=dnw, min_size=5)
#'
#' #view the results
#' netGSAresults(dnw)
#'
#'
#' @importFrom stats p.adjust
#' @import igraph
#' @importFrom netgsa NetGSA
#' @export
runNetGSA <- function(object,
                      min_size=5){

  ##test for proper input
  if(!inherits(object, "DNEAobj")) stop('the input object should be of class "DNEAobj"!')
  if(min_size <1) stop("min_size parameter should be a positive integer greater than zero!")

  ##set input variables
  adjacency_matrices <- list(list(adjacencyMatrix(x=object, weighted=TRUE)[[1]]),
                             list(adjacencyMatrix(x=object, weighted=TRUE)[[2]]))

  expression_data <- expressionData(x=object, assay="log_input_data")
  data_groups <- ifelse(networkGroupIDs(object) == networkGroups(object)[1], 1, 2)
  subnetworks <- as.matrix(subnetworkMembership(object))

  ##filter subnetworks to only include those greater than or equal to min_size
  sub_membership <- subnetworkMembership(object)
  filtered_subnetworks <- as.matrix(sub_membership[rowSums(sub_membership) >= min_size, ])

  ##run netgsa
  netgsa_results <- NetGSA(A=adjacency_matrices,
                           x=expression_data,
                           group=data_groups,
                           pathways=filtered_subnetworks,
                           lklMethod="REML",
                           minsize=min_size)

  ##add netGSA results to Node list
  nodeList(object)[["mean1"]] <- as.vector(netgsa_results$beta[[1]])
  nodeList(object)[["mean2"]] <- as.vector(netgsa_results$beta[[2]])
  nodeList(object)[["meanchange"]] <- netgsa_results$beta[[2]] - netgsa_results$beta[[1]]
  nodeList(object)[["mc.notes"]] <- paste(networkGroups(object)[[2]], 'over', networkGroups(object)[[1]])

  ##concatenate netGSA summary output
  res <- data.frame(CCsummary(object)[CCsummary(object)$number_of_nodes >= min_size, ])
  res <- res[!grepl("independent", res$subnetworks),]
  res[["NetGSA_pval"]] <- netgsa_results$results$pval
  res[["NetGSA_pFDR"]] <- netgsa_results$results$pFdr
  res <- res[order(res$NetGSA_pFDR),]

  ##rename subnetworks
  cluster_names <- CCsummary(object)[-match("independent", CCsummary(object)$subnetworks), ]
  new_cluster_order <- data.frame(subnetworks=c(res$subnetworks, cluster_names$subnetworks[-match(res$subnetworks, cluster_names$subnetworks)]))

  ##update consensus subnetworks
  nodeList(object)$membership <- match(nodeList(object)$membership, as.numeric(gsub('subnetwork','',new_cluster_order$subnetworks)))
  nodeList(object)$membership[is.na(nodeList(object)$membership)] <- "independent"

  ##update res subnetwork names
  res$subnetworks <- paste0("subnetwork", seq(1, nrow(res)))

  ##update cluster summary subnetwork names
  cluster_names <- CCsummary(object)
  rownames(cluster_names) <- cluster_names$subnetworks
  cluster_names <- cluster_names[new_cluster_order$subnetworks, ]
  rownames(cluster_names) <- NULL
  cluster_names$subnetworks <- paste0("subnetwork", seq(1, nrow(cluster_names)))
  cluster_names <- rbind(cluster_names, CCsummary(object)[match("independent", CCsummary(object)$subnetworks), ])

  ##update subnetwork_membership matrix
  sub_membership <- sub_membership[new_cluster_order$subnetworks, ]
  rownames(sub_membership) <- paste0("subnetwork", seq(1, nrow(sub_membership)))

  ##update DNEAobject
  netGSAresults(object) <- res
  CCsummary(object) <- cluster_names
  subnetworkMembership(object) <- sub_membership

  ##check valid object
  validObject(object)
  return(object)
}
