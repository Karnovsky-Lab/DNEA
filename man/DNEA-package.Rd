% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/start-here.R
\docType{package}
\name{DNEA-package}
\alias{DNEA-package}
\title{DNEA package for joint estimation of biological networks}
\description{
DNEA is designed to take as input non-normalized, non-transformed -omics data
(ie. metabolomics, lipidomics, proteomics) and normalize the expression values, optimize the regularization
parameter, and jointly estimate the biological networks across a user-specified experimental
condition. The resulting networks may then be clustered to identify metabolic modules,
aka subnetworks, within the larger networks. These subnetworks are then tested for enrichment across the experimental
condition using the \code{\link[netgsa:NetGSA]{netgsa::NetGSA()}} algorithm.
}
\details{
The main workflow contains the following functions:
\enumerate{
\item \strong{\code{\link{createDNEAobject}}}
\item \strong{\code{\link{BICtune}}}
\item \strong{\code{\link{stabilitySelection}}}
\item \strong{\code{\link{getNetworks}}}
\item \strong{\code{\link{clusterNet}}}
\item \strong{\code{\link{runNetGSA}}}}

A more descriptive workflow can be viewed in the package vignette. This can be accessed by running \code{vignette("DNEA")} in the console. \cr

The source code is available at the Karnovsky lab Github page: \url{https://github.com/Karnovsky-Lab/DNEAdev}
}
\keyword{package}