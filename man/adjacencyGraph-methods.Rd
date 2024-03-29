% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-generics.R, R/all-methods.R
\name{adjacencyGraph}
\alias{adjacencyGraph}
\alias{adjacencyGraph,DNEAobj-method}
\alias{adjacencyGraph,consensusClusteringResults-method}
\title{Retrieve the adjacency graph for the case, control, or joint network}
\usage{
adjacencyGraph(x, graph)

\S4method{adjacencyGraph}{DNEAobj}(x, graph)

\S4method{adjacencyGraph}{consensusClusteringResults}(x, graph)
}
\arguments{
\item{x}{A \code{\link{DNEAobj}} or \code{consensusClusteringResults} object}

\item{graph}{A character string indicating which of the adjacency
graphs to return. Values can be "joint_graph" for the whole graph object,
or one of the group values returned by \code{\link{networkGroups}}}
}
\value{
An \code{\link{igraph}} graph object corresponding to
the specified adjacency graph.
}
\description{
The function  returns the adjacency graph made for the case,
control, or joint network determined via \code{\link{clusterNet}}.
}
\examples{
#import example data
data(dnw)

adjacencyGraph(dnw, graph = "DM:case")
}
\seealso{
\code{\link{clusterNet}}
}
\author{
Christopher Patsalis
}
