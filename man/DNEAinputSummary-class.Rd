% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-classes.R, R/all-methods.R
\docType{class}
\name{DNEAinputSummary-class}
\alias{DNEAinputSummary-class}
\alias{DNEAinputSummary}
\alias{show,DNEAinputSummary-method}
\title{DNEAinputSummary}
\usage{
\S4method{show}{DNEAinputSummary}(object)
}
\arguments{
\item{object}{A DNEAinputSummary object}
}
\value{
A DNEAinputSummary object

A summary of the input data to \code{\link{createDNEAobject}}.
}
\description{
An s4 class to represent the results from diagnostic testing
on the input data
}
\section{Functions}{
\itemize{
\item \code{show(DNEAinputSummary)}: This function will display the number of samples, number of features,
and diagnostics values of the input dataset to a \code{\link{DNEAobj}} object.

}}
\section{Slots}{

\describe{
\item{\code{num_samples}}{a single-value numeric vector corresponding to the
number of samples in the dataset}

\item{\code{num_features}}{a single-value numeric vector corresponding to the
number of features in the dataset}

\item{\code{diagnostic_values}}{a 3x3 data frame with the diagnostic values
calculated via \code{\link{createDNEAobject}}}
}}

\seealso{
\code{\link{createDNEAobject}}

\code{\link{createDNEAobject}},
\code{\link{aggregateFeatures}}
}
\author{
Christopher Patsalis
}
