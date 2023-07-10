#################################################################################################################################################
#* DNEA package organization *#                                                                                                                 #
###############################                                                                                                                 #
#                                                                                                                                               #
# start-here.R.....................This file - contains the package legend as well as package and data documentation                            #
# createDNEAobject.R...............contains the exported createDNEAobject() function that initiates the DNEA workflow. The internal helper      #
#                                  functions called by createDNEAobject (ie. restructure_input_data(), dataDiagnostics(), metabDE()) are also   #
#                                  located in this file                                                                                         #
# collapseNodes-fun.R..............contains the exported wrapper function(reduceFeatures()) for node collapsing as well as the internal         #
#                                  helper functions called by said function (ie. collapseNodes_cor(), collapseNodes_knowledge(),                #
#                                  and collapseNodes_hybrid())                                                                                  #
# MainFunctions.R..................contains the exported wrapper functions for each step of DNEA (ie. BICtune(), stabilitySelection(),          #
#                                  getNetworks(), clusterNet(), runNetGSA())                                                                    #
# utilities-internal.R.............contains internal utility functions used by exported functions (ie. split_by_condition)                      #
# utilities-external.R.............contains exported utility functions that adds functionality for the user but not part of                     #
#                                  the DNEA workflow (ie. includeMetadata(), plotNetworks(), getNetworkFiles(), filterNetworks())               #
# all-classes.R....................contains the custom s4 classes and validator functions                                                       #
# all-generics.R...................contains all s4 generics                                                                                     #
# all-methods.R....................contains all s4 custom methods                                                                               #
# JSEM.R...........................contains all of the internal functions for joint estimation of biological networks (ie. functions            #
#                                  called by BICtune(), stabilitySelection(), and getNetworks())                                                #
# preprocess_lib.R.................contains the internal functions necessary for consensus clustering                                           #
#################################################################################################################################################

#' DNEA packge for joint estimation of biological networks
#'
#' DNEA is designed to take as input non-normalized, non-transformed expression data -omics data
#' (metabolomics, lipidomics, protiomics) and normalize the data, optimize the regularization
#' parameter and fit jointly estimate the biological networks across a specified experimental
#' condition. After Network construction, the user may cluster them using to identify metabolic modules,
#' or subnetworks, within the larger networks and subsequently test for enrichment across the experimental
#' condition of said subnetworks using the \code{\link{netgsa, NetGSA}} algorithm.
#'
#' The main workflow contains the following functions:
#' \enumerate{
#' \item \strong{\code{\link{createDNEAobject}}}
#' \item \strong{\code{\link{BICtune}}}
#' \item \strong{\code{\link{stabilitySelection}}}
#' \item \strong{\code{\link{getNetworks}}}
#' \item \strong{\code{\link{clusterNet}}}
#' \item \strong{\code{\link{runNetGSA}}}}
#'
#' A more descriptive workflow can be view via the package vignette. This can be accessed by running \code{vignette("DNEA")}. \cr
#'
#' The source code is available at the Karnovsky lab Github page: \url{https://github.com/Karnovsky-Lab/DNEAdev}
#'
#' @docType package
#' @name DNEA-package
#' @keywords package

#' "TEDDY" example data from xxx
#'
#' This data is an \emph{n x m} expression matrix corresponding to a curated list of metabolites from "The Environmental Determinants of
#' Diabetes in the Young" clinical trial.
#'
#' @format ## `TEDDY`
#' A data frame with 100 rows and 144 columns. Each row corresponds to a sample, and each column corresponds to a metabolite.
#'
#' @references Lee HS, Burkhardt BR, McLeod W, Smith S, Eberhard C, Lynch K, Hadley D, Rewers M, Simell O, She JX, Hagopian B, Lernmark A, Akolkar B, Ziegler AG, Krischer JP; TEDDY study group. Biomarker discovery study design for type 1 diabetes in The Environmental Determinants of Diabetes in the Young (TEDDY) study. Diabetes Metab Res Rev. 2014 Jul;30(5):424-34. doi: 10.1002/dmrr.2510. PMID: 24339168; PMCID: PMC4058423. \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4058423/}
#' @source The raw data can be downloaded from the Metabolomics workbench under study ID \strong{ST001386}:
#' \url{https://www.metabolomicsworkbench.org/data/DRCCStudySummary.php?Mode=SetupRawDataDownload&StudyID=ST001386}
"TEDDY"
