################################################################################
#* DNEA package organization *#                                                #                                                                #
###############################                                                #                                                                                         #
# start-here.R.....................This file - contains the package legend as  #
#                                  well as package and data documentation      #
# initiator.R......................contains the exported createDNEAobject()    #
#                                  function that initiates the DNEA workflow.  #
#                                  The internal helper functions called by     #
#                                  createDNEAobject (ie.                       #
#                                  restructure_input_data(), dataDiagnostics(),#
#                                  metabDE()) are also located in this file    #                                                                                    #
# reduce-features.R................contains the exported wrapper function      #
#                                  aggregateFeatures() for node aggregating as #
#                                  well as the internal helper functions called#
#                                  by said function (ie. collapseNodes_cor(),  #
#                                  collapseNodes_knowledge(), and              #
#                                  collapseNodes_hybrid())                     #
# primary.R........................contains the exported wrapper functions for #
#                                  each step of DNEA (ie. BICtune(),           #
#                                  stabilitySelection(), getNetworks(),        #
#                                  clusterNet(), runNetGSA())                  #
# utilities-internals.R............contains internal utility functions used by #
#                                  exported functions (ie. split_by_condition) #
# utilities-exported.R.............contains exported utility functions that    #
#                                  adds functionality for the user but not part#
#                                  of the DNEA workflow (ie. includeMetadata(),#
#                                  plotNetworks(), getNetworkFiles(),          #
#                                  filterNetworks())                           #
# all-classes.R....................contains the custom s4 classes and validator#
#                                  functions                                   #
# all-generics.R...................contains all s4 generics                    #                                                                   #
# all-methods.R....................contains all s4 custom methods              #                                                                #
# JSEM-internals.R.................contains all of the internal functions for  #
#                                  joint estimation of biological networks (ie.#
#                                  functions called by BICtune(),              #
#                                  stabilitySelection(), and getNetworks())    #
# clustering-internals.R...........contains the internal functions necessary   #
#                                  for consensus clustering                    #
################################################################################
#' @docType package
#' @name DNEA-package
#' @rdname DNEA-package
#' @aliases DNEA-package DNEA
#' @section Primary Components:
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
#' A more descriptive workflow can be viewed in the package vignette.
#' This can be accessed by running \code{vignette("DNEA")} in the console.
#'
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL


#' Example expresion data set from The Environmental Determinants
#' of Diabetes in the Young (TEDDY) clinical trial
#'
#' This data is an \emph{m x n} expression matrix corresponding to a
#' curated list of metabolites from "The Environmental Determinants of
#' Diabetes in the Young" clinical trial. The data was downloaded from
#'
#'
#' @returns An \emph{m x n} expression matrix of metabolomics data from the
#' TEDDY dataset
#' @docType data
#' @keywords datasets
#' @name TEDDY
#' @usage data("TEDDY")
#' @format A data frame with 134 rows and 322 columns. Each row corresponds
#' to a unique metabolite, and each column corresponds to a sample
#' @references Lee HS, Burkhardt BR, McLeod W, Smith S, Eberhard C, Lynch K,
#' Hadley D, Rewers M, Simell O, She JX, Hagopian B, Lernmark A, Akolkar B,
#' Ziegler AG, Krischer JP; TEDDY study group. Biomarker discovery study
#' design for type 1 diabetes in The Environmental Determinants of Diabetes
#' in the Young (TEDDY) study. Diabetes Metab Res Rev. 2014 Jul;30(5):424-34.
#' doi: 10.1002/dmrr.2510. PMID: 24339168;
#' PMCID:
#' \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4058423/}{PMC4058423}
#'
#' @source This data is available at the NIH Common Fund's National
#' Metabolomics Data Repository (NMDR) website,
#' \href{https://www.metabolomicsworkbench.org}{the Metabolomics Workbench},
#'  where it has been assigned
#' Project ID
#' \href{https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR000950}{PR000950}
#'  and Study ID ST001386. The data can be accessed
#' directly via it's Project DOI: \url{10.21228/M8WM4P}. This work is
#' supported by NIH grant, U2C- DK119886.
"TEDDY"

#' Sample meta data for the The Environmental Determinants
#' of Diabetes in the Young (TEDDY) clinical trial
#'
#' This is a data frame containing metadata for the samples
#' in the corresponding \code{\link{TEDDY}} example data from
#' "The Environmental Determinants of Diabetes in the Young"
#' clinical trial.
#'
#' @returns A data frame containing the sample metadata for
#' the TEDDY metabolomics study
#' @docType data
#' @keywords datasets
#' @name T1Dmeta
#' @usage data("T1Dmeta")
#' @format A data frame with 322 rows and 7 columns. Each row corresponds
#' to a sample, and each column corresponds to:
#' \describe{
#' \item{subject}{The individual patient}
#' \item{Endpoint1}{The age of the case subject in days when this
#' sample was collected}
#' \item{Endpoint2}{The age of the control subject in days when
#' this sample was collected}
#' \item{Age}{The age of the subject in days when this sample
#' was collected}
#' \item{Sex}{The sex of the subject}
#' \item{sample}{The name of this sample}
#' \item{group}{A variable indicating whether or not this sample
#' is part of the T1D case or T1D control group}}
#'
#' @references Lee HS, Burkhardt BR, McLeod W, Smith S, Eberhard C,
#' Lynch K, Hadley D, Rewers M, Simell O, She JX, Hagopian B,
#' Lernmark A, Akolkar B, Ziegler AG, Krischer JP; TEDDY study group.
#' Biomarker discovery study design for type 1 diabetes in The
#' Environmental Determinants of Diabetes in the Young (TEDDY) study.
#' Diabetes Metab Res Rev. 2014 Jul;30(5):424-34. doi: 10.1002/dmrr.2510.
#' PMID: 24339168; PMCID: PMC4058423.
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4058423/}
#'
#' @source The raw data can be downloaded from the Metabolomics workbench
#' under study ID \strong{ST001386}:
#' \url{https://www.metabolomicsworkbench.org/data/DRCCStudySummary.php?Mode=SetupRawDataDownload&StudyID=ST001386}
"T1Dmeta"

#' Feature meta data for the The Environmental Determinants
#' of Diabetes in the Young (TEDDY) clinical trial
#'
#' This is a data frame containing metadata for the metabolites
#' in the corresponding \code{\link{TEDDY}} example data from
#' "The Environmental Determinants of Diabetes in the Young"
#' clinical trial.
#'
#' @returns A data frame containing the metabolite metadata
#' for the TEDDY metabolomics study
#' @docType data
#' @keywords datasets
#' @name metab_data
#' @usage data("metab_data")
#' @format A data frame with 134 rows and 3 columns. Each row corresponds
#' to a metabolite, and each column corresponds to:
#' \describe{
#' \item{variable_id}{The metabolite name}
#' \item{mz}{The mass/charge ratio for a given metabolite}
#' \item{rt}{The retention time for a given metabolite}}
#'
#' @references Lee HS, Burkhardt BR, McLeod W, Smith S, Eberhard C,
#' Lynch K, Hadley D, Rewers M, Simell O, She JX, Hagopian B,
#' Lernmark A, Akolkar B, Ziegler AG, Krischer JP; TEDDY study group.
#' Biomarker discovery study design for type 1 diabetes in The
#' Environmental Determinants of Diabetes in the Young (TEDDY) study.
#' Diabetes Metab Res Rev. 2014 Jul;30(5):424-34. doi: 10.1002/dmrr.2510.
#' PMID: 24339168; PMCID: PMC4058423.
#' \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4058423/}
#'
#' @source The raw data can be downloaded from the Metabolomics workbench
#' under study ID \strong{ST001386}:
#' \url{https://www.metabolomicsworkbench.org/data/DRCCStudySummary.php?Mode=SetupRawDataDownload&StudyID=ST001386}
"metab_data"

#' Example results for DNEA
#'
#' "dnw" is a DNEA object containing the results for the full
#' DNEA workflow on the \code{\link{TEDDY}} example data. The
#' exact workflow to produce these results can be replicated by
#' following the package vignette accessed by entering
#' browseVignettes("DNEA") in the console. 1000 replicates were
#' performed during stability selection \emph{with} the
#' subsampling protocol. The lambda value used during joint
#' estimation was aproximated as
#' \deqn{\lambda = \sqrt{ \ln (num. features) / num. samples}}{ lambda = sqrt(ln(num. features) / num. samples)}
#'
#'
#' @returns A \code{\link[=DNEA-class]{DNEA}} object
#' containing the results of a DNEA experiment.
#' @docType data
#' @keywords datasets
#' @name dnw
#' @usage data("dnw")
#' @format A DNEA results object after completing a DNEA experiment.
#' @source The data the results of the full DNEA workflow performed using
#' the \code{\link{TEDDY}} example data, as described above.
"dnw"


