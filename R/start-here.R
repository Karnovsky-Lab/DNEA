#################################################################################################################################################
#* DNEAdev package organization *#                                                                                                                 #
###############################                                                                                                                 #
#                                                                                                                                               #
# start-here.R.....................This file - contains the package legend as well as package and data documentation                            #
# createDNEAobject.R...............contains the exported createDNEAobject() function that initiates the DNEAdev workflow. The internal helper      #
#                                  functions called by createDNEAobject (ie. restructure_input_data(), dataDiagnostics(), metabDE()) are also   #
#                                  located in this file                                                                                         #
# collapseNodes-fun.R..............contains the exported wrapper function(reduceFeatures()) for node collapsing as well as the internal         #
#                                  helper functions called by said function (ie. collapseNodes_cor(), collapseNodes_knowledge(),                #
#                                  and collapseNodes_hybrid())                                                                                  #
# MainFunctions.R..................contains the exported wrapper functions for each step of DNEAdev (ie. BICtune(), stabilitySelection(),          #
#                                  getNetworks(), clusterNet(), runNetGSA())                                                                    #
# utilities-internal.R.............contains internal utility functions used by exported functions (ie. split_by_condition)                      #
# utilities-external.R.............contains exported utility functions that adds functionality for the user but not part of                     #
#                                  the DNEAdev workflow (ie. includeMetadata(), plotNetworks(), getNetworkFiles(), filterNetworks())               #
# all-classes.R....................contains the custom s4 classes and validator functions                                                       #
# all-generics.R...................contains all s4 generics                                                                                     #
# all-methods.R....................contains all s4 custom methods                                                                               #
# JSEM.R...........................contains all of the internal functions for joint estimation of biological networks (ie. functions            #
#                                  called by BICtune(), stabilitySelection(), and getNetworks())                                                #
# preprocess_lib.R.................contains the internal functions necessary for consensus clustering                                           #
#################################################################################################################################################

#' The example data from xxx
#'
#' This data is an \emph{m x n} expression matrix corresponding to a curated list of metabolites from "The Environmental Determinants of
#' Diabetes in the Young" clinical trial.
#'
#' @returns An \emph{m x n} expression matrix of metabolomics data from the TEDDY dataset
#' @docType data
#' @keywords datasets
#' @name TEDDY
#' @usage data("TEDDY")
#' @format A data frame with 100 rows and 144 columns. Each row corresponds to a sample, and each column
#'         corresponds to a metabolite.
#' @references Lee HS, Burkhardt BR, McLeod W, Smith S, Eberhard C, Lynch K, Hadley D, Rewers M, Simell O, She JX, Hagopian B, Lernmark A, Akolkar B, Ziegler AG, Krischer JP; TEDDY study group. Biomarker discovery study design for type 1 diabetes in The Environmental Determinants of Diabetes in the Young (TEDDY) study. Diabetes Metab Res Rev. 2014 Jul;30(5):424-34. doi: 10.1002/dmrr.2510. PMID: 24339168; PMCID: PMC4058423. \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4058423/}
#' @source This data is available at the NIH Common Fund's National Metabolomics Data Repository (NMDR) website,
#'         the Metabolomics Workbench, \url{https://www.metabolomicsworkbench.org}, where it has been assigned
#'         Project ID PR000950 and Study ID ST001386. The data can be accessed directly via it's
#'         Project DOI: \url{10.21228/M8WM4P}. This work is supported by NIH grant, U2C- DK119886.
"TEDDY"

#' The meta data file for the TEDDY metabolomics data
#'
#' This data is a data frame containing metadata for the samples in the corresponding \code{\link{TEDDY}} example data from "The Environmental Determinants of
#' Diabetes in the Young" clinical trial.
#'
#' @returns A data frame containing the metadata for the TEDDY metabolomics study
#' @docType data
#' @keywords datasets
#' @name T1Dmeta
#' @usage data("T1Dmeta")
#' @format A data frame with 322 rows and 9 columns. Each row corresponds to a sample, and each column corresponds to:
#' \describe{
#'   \item{subject}{The individual patient}
#'   \item{Endpoint1}{The age of the case subject in days when this sample was collected}
#'   \item{Endpoint2}{The age of the control subject in days when this sample was collected}
#'   \item{Age}{The age of the subject in days when this sample was collected}
#'   \item{Sex}{The sex of the subject}
#'   \item{sample}{The name of this sample}
#'   \item{group}{A variable indicating whether or not this sample is part of the T1D case or T1D control group}
#' }
#' @references Lee HS, Burkhardt BR, McLeod W, Smith S, Eberhard C, Lynch K, Hadley D, Rewers M, Simell O, She JX, Hagopian B, Lernmark A, Akolkar B, Ziegler AG, Krischer JP; TEDDY study group. Biomarker discovery study design for type 1 diabetes in The Environmental Determinants of Diabetes in the Young (TEDDY) study. Diabetes Metab Res Rev. 2014 Jul;30(5):424-34. doi: 10.1002/dmrr.2510. PMID: 24339168; PMCID: PMC4058423. \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4058423/}
#' @source The raw data can be downloaded from the Metabolomics workbench under study ID \strong{ST001386}:
#' \url{https://www.metabolomicsworkbench.org/data/DRCCStudySummary.php?Mode=SetupRawDataDownload&StudyID=ST001386}
"T1Dmeta"

#' Example results for DNEAdev
#'
#' This data is an example DNEAresults object after performing DNEAdev. The experiment was performed using the \code{\link{TEDDY}}
#' data with 4 reps of stability selection, no subsampling, and default parameters everywhere else. The seed was set to 417.
#'
#' @returns A \code{DNEAresults} object containing the results of a DNEAdev experiment
#' @docType data
#' @keywords datasets
#' @name dnw
#' @usage data("dnw")
#' @format A DNEAdev results object after completing a DNEAdev experiment.
#' @source The data is stored in the \code{\link{DNEAdev}} package and can be accessed by running data(dnw) in the console.
"dnw"

#' DNEAdev package for joint estimation of biological networks
#'
#' DNEAdev is designed to take as input non-normalized, non-transformed -omics data
#' (ie. metabolomics, lipidomics, proteomics) and normalize the expression values, optimize the regularization
#' parameter, and jointly estimate the biological networks across a user-specified experimental
#' condition. The resulting networks may then be clustered to identify metabolic modules,
#' aka subnetworks, within the larger networks. These subnetworks are then tested for enrichment across the experimental
#' condition using the \code{\link[netgsa:NetGSA]{netgsa::NetGSA()}} algorithm.
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
#' A more descriptive workflow can be viewed in the package vignette. This can be accessed by running \code{vignette("DNEAdev")} in the console. \cr
#'
#' The source code is available at the Karnovsky lab Github page: \url{https://github.com/Karnovsky-Lab/DNEAdev}
#'
#' @docType package
#' @name DNEAdev
#' @aliases DNEAdev
#' @keywords package
NULL
