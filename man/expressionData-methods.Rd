% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-generics.R, R/all-methods.R
\name{expressionData}
\alias{expressionData}
\alias{expressionData,DNEAobj-method}
\title{Access expression data within a DNEAobj object,}
\usage{
expressionData(x, normalized)

\S4method{expressionData}{DNEAobj}(x, normalized = TRUE)
}
\arguments{
\item{x}{A \code{\link{DNEAobj}} object}

\item{normalized}{A boolean indicating whether the normalized or
original input data should be returned}
}
\value{
The expression matrix specified by the user.
}
\description{
This function accesses the expression data stored in the @assays
slot of the \code{\link{DNEAobj}} object. The output is an \emph{n x m} matrix
with one row for each sample and one column for each feature in the data.
}
\examples{
#import example data
data(dnw)

expressionData(dnw, normalized = TRUE)
}
\seealso{
\code{\link{createDNEAobject}},\code{\link{aggregateFeatures}}
}
\author{
Christopher Patsalis
}