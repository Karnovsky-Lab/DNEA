% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-generics.R, R/all-methods.R
\name{adjacencyMatrix}
\alias{adjacencyMatrix}
\alias{adjacencyMatrix,DNEAobj-method}
\title{Retrieve the weighted or unweighted adjacency matrix}
\usage{
adjacencyMatrix(x, weighted)

\S4method{adjacencyMatrix}{DNEAobj}(x, weighted = FALSE)
}
\arguments{
\item{x}{A \code{\link{DNEAobj}} object}

\item{weighted}{A boolean indicating whether or not
to select the weighted or thresholded (unweighted)
adjacency matrix}
}
\value{
A matrix corresponding to the adjacency matrix specified.
}
\description{
The function takes as input a \code{\link{DNEAobj}} object and returns the
weighted or unweighted adjacency matrix determined via
\code{\link{getNetworks}}.
}
\examples{
#import example data
data(dnw)

adjacencyMatrix(dnw, weighted = TRUE)
}
\seealso{
\code{\link{getNetworks}}
}
\author{
Christopher Patsalis
}