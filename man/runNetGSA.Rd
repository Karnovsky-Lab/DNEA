% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/primary.R
\name{runNetGSA}
\alias{runNetGSA}
\title{Identify metabolic modules that are enriched across
experimental conditions using NetGSA}
\usage{
runNetGSA(object, min_size = 5, assay = "log_input_data")
}
\arguments{
\item{object}{A \code{\link{DNEAobj}}.}

\item{min_size}{The minimum size of a given metabolic
module for to be tested for enrichment across the
experimental condition.}

\item{assay}{A character string indicating which expression assay to
use for analysis. The default is the "log_input-data" assay that is
created during \code{\link{createDNEAobject}}.}
}
\value{
A \code{\link{DNEAobj}} object after populating the @netGSA
slot. A summary of the NetGSA results can be viewed
using \code{\link{netGSAresults}}.
}
\description{
This function performs pathway enrichment analysis on the metabolic
modules identified via \code{\link{clusterNet}} using the
\code{\link[netgsa:NetGSA]{netgsa::NetGSA()}} algorithm.
}
\examples{
#dnw is a DNEAobj with the results generated for the example data
#accessed by running data(TEDDY) in the console. The workflow
#for this data can be found in the vignette accessed by
#running browseVignettes("DNEA") in the console.
data(dnw)

#perform pathway enrichment analysis using netGSA
dnw <- runNetGSA(object=dnw, min_size=5)

#view the results
netGSAresults(dnw)


}
\references{
Hellstern M, Ma J, Yue K, Shojaie A.
netgsa: Fast computation and interactive visualization for
topology-based pathway enrichment analysis.
PLoS Comput Biol. 2021 Jun 11;17(6):e1008979.
doi: 10.1371/journal.pcbi.1008979. PMID: 34115744;
PMCID: PMC8221786
\url{https://pubmed.ncbi.nlm.nih.gov/34115744/}
}
\seealso{
\code{\link{netGSAresults}}
\code{\link{clusterNet}}
\code{\link[netgsa:NetGSA]{NetGSA}}
}
\author{
Christopher Patsalis
}
