% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-generics.R, R/all-methods.R
\name{CCsummary}
\alias{CCsummary}
\alias{CCsummary,DNEA-method}
\title{Retrieves the summary results of consensus clustering}
\usage{
CCsummary(x)

\S4method{CCsummary}{DNEA}(x)
}
\arguments{
\item{x}{A \code{\link{DNEA}} object.}
}
\value{
A data frame summary of the consensus clustering
results from DNEA.
}
\description{
The function takes as input a \code{\link{DNEA}} object and
returns a summary  of the results of consensus clustering
stored in the consensus_clustering slot as a
\code{\link{consensusClusteringResults}} object.
}
\examples{
#dnw is a \code{\link{DNEA}} object with the results
#generated for the example data accessed by running
#data(TEDDY) in the console. The workflow for this data
#can be found in the vignette accessed by running
#browseVignettes("DNEA") in the console.
data("dnw")

CCsummary(dnw)
}
\seealso{
\code{\link{clusterNet}}
}
\author{
Christopher Patsalis
}
