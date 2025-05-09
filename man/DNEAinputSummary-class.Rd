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
on the input data to a \code{\link[=DNEA-class]{DNEA}}.
}
\section{Functions}{
\itemize{
\item \code{show(DNEAinputSummary)}: This function will display the number of samples, number of features,
and diagnostics values of the input data set to a
\code{\link{DNEA}} object.

}}
\section{Slots}{

\describe{
\item{\code{num_samples}}{a single-value numeric vector corresponding
to the number of samples in the data set.}

\item{\code{num_features}}{a single-value numeric vector corresponding
to the number of features in the data set}

\item{\code{diagnostic_values}}{a 3x3 data frame with the diagnostic
values calculated via \code{\link{createDNEAobject}}.}
}}

\examples{
#load example data
data(TEDDY)
data(T1Dmeta)

#make sure metadata and expression data are in same order
T1Dmeta <- T1Dmeta[colnames(TEDDY),]

#create group labels
group_labels <- T1Dmeta$group
names(group_labels) <- rownames(T1Dmeta)

#initiate DNEA object
dnw <- createDNEAobject(project_name = "test", expression_data = TEDDY,
                            group_labels = group_labels)

datasetSummary(dnw)
}
\seealso{
\code{\link{createDNEAobject}}

\code{\link{createDNEAobject}},
\code{\link{aggregateFeatures}}
}
\author{
Christopher Patsalis
}
