% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all-generics.R, R/all-methods.R
\name{sampleNames}
\alias{sampleNames}
\alias{sampleNames,DNEA-method}
\title{Retrieve the sample names from the metadata slot.}
\usage{
sampleNames(x, original = FALSE)

\S4method{sampleNames}{DNEA}(x, original = FALSE)
}
\arguments{
\item{x}{A \code{\link{DNEA}} object.}

\item{original}{"TRUE" returns the original sample names
and "FALSE" returns the sample names that have been
modified to avoid errors as a result of special characters
using \code{\link[janitor:make_clean_names]{make_clean_names}}.}
}
\value{
A character vector of sample names.
}
\description{
This function accesses the sample names stored in the
metadata slot of the \code{\link{DNEA}} object.
}
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

sampleNames(dnw)
}
\seealso{
\code{\link{createDNEAobject}}
}
\author{
Christopher Patsalis
}
