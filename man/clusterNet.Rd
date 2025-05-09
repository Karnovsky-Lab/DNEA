% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/primary.R
\name{clusterNet}
\alias{clusterNet}
\title{Identify metabolic modules within the biological networks using a
consensus clustering approach}
\usage{
clusterNet(object, tau = 0.5, max_iterations = 5, verbose = TRUE)
}
\arguments{
\item{object}{A \code{\link[=DNEA-class]{DNEA}} object.}

\item{tau}{The \% agreement among the clustering algorithms
for a node to be included in a sub network.}

\item{max_iterations}{The maximum number of replicates of
consensus clustering to be performed if consensus is not
reached.}

\item{verbose}{TRUE/FALSE whether a progress bar should be
displayed in the console.}
}
\value{
A \code{\link[=DNEA-class]{DNEA}} object containing sub network
determinations for the nodes within the input network. A summary of the
consensus clustering results can be viewed using \code{\link{CCsummary}}.
Sub network membership for each node can be found in the "membership"
column of the node list, which can be accessed using \code{\link{nodeList}}.
}
\description{
This function clusters the jointly estimated adjacency matrices
constructed using \code{\link{getNetworks}} via the consensus clustering
approach described in Ma et al. (\emph{Please see the
\strong{\emph{Details}} section for more information}) to identify
metabolic modules, aka sub networks, present in the larger networks.
Only sub networks with consensus that meets or exceeds tau are
identified as real.
}
\details{
Seven clustering algorithms from the \code{\link{igraph}} package
are utilized in this consensus clustering approach:
\enumerate{
\item \code{\link[igraph:cluster_edge_betweenness]{cluster_edge_betweenness}}
\item \code{\link[igraph:cluster_fast_greedy]{cluster_fast_greedy}}
\item \code{\link[igraph:cluster_infomap]{cluster_infomap}}
\item \code{\link[igraph:cluster_label_prop]{cluster_label_prop}}
\item \code{\link[igraph:cluster_louvain]{cluster_louvain}}
\item \code{\link[igraph:cluster_walktrap]{cluster_walktrap}}
\item \code{\link[igraph:cluster_leading_eigen]{cluster_leading_eigen}}}

For each iteration, node membership in a respective cluster is
compared across the algorithms, and only the nodes with tau \%
agreement for a given cluster are kept. A new adjacency graph is
then created and clustering is performed again. This occurs iteratively
until consensus on is reached stable sub networks or the specified
"max_iterations" is reached
\emph{(Please see references for more details)}.
}
\examples{
#dnw is a \code{\link[=DNEA-class]{DNEA}} object with the results
#generated for the example data accessed by running
#data(TEDDY) in the console. The workflow for this data
#can be found in the vignette accessed by running
#browseVignettes("DNEA") in the console.
data(dnw)

#identify metabolic modules via consensus clustering
dnw <- clusterNet(object=dnw, tau=0.5, max_iterations=5)

#we can also plot the subnetworks
plotNetworks(object=dnw, type="sub_networks", subtype=1)

}
\references{
Ma J, Karnovsky A, Afshinnia F, Wigginton J, Rader DJ, Natarajan L,
Sharma K, Porter AC, Rahman M, He J, Hamm L, Shafi T, Gipson D,
Gadegbeku C, Feldman H, Michailidis G, Pennathur S.
Differential network enrichment analysis reveals novel lipid pathways
in chronic kidney disease. Bioinformatics. 2019 Sep 15;35(18):3441-3452.
doi: 10.1093/bioinformatics/btz114. PMID: 30887029; PMCID: PMC6748777.
\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6748777/}
}
\seealso{
\code{\link{plotNetworks}}
}
\author{
Christopher Patsalis
}
