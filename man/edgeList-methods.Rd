% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-generics.R, R/all-methods.R
\name{edgeList}
\alias{edgeList}
\alias{edgeList<-}
\alias{edgeList,DNEA-method}
\alias{edgeList<-,DNEA-method}
\title{Access the edge list}
\usage{
edgeList(x)

edgeList(x) <- value

\S4method{edgeList}{DNEA}(x)

\S4method{edgeList}{DNEA}(x) <- value
}
\arguments{
\item{x}{a \code{\link{DNEA}} object.}

\item{value}{a data frame of edges in the network.}
}
\value{
A data frame corresponding to the edge list
determined by DNEA.
}
\description{
The function takes as input a \code{\link{DNEA}} object and
returns the edge list created by the \code{\link{getNetworks}}
function.
}
\examples{
#dnw is a \code{\link{DNEA}} object with the results
#generated for the example data accessed by running
#data(TEDDY) in the console. The workflow for this data
#can be found in the vignette accessed by running
#browseVignettes("DNEA") in the console.
data("dnw")

edgeList(dnw)
}
\seealso{
\code{\link{getNetworks}},\code{\link{filterNetworks}},
\code{\link{getNetworkFiles}}
}
\author{
Christopher Patsalis
}
