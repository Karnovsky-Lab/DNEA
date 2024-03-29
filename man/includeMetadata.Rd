% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities-exported.R
\name{includeMetadata}
\alias{includeMetadata}
\title{Add additional metadata to the DNEAobj object}
\usage{
includeMetadata(object, type = c("sample", "feature"), metadata)
}
\arguments{
\item{object}{A \code{\link{DNEAobj}} object}

\item{type}{sample or feature metadata}

\item{metadata}{a data.frame containing metadata to add}
}
\value{
A DNEAobj object with the specified additions
}
\description{
This function will take additional metadata and add it to the specified
dataframe in the metadata slot. \strong{\emph{NOTE:}} The rownames of the
new metadata must match the order of the input sample names or feature
names, respectively
}
\examples{
#import example data
data(dnw)
data(T1Dmeta)

#make sure metadata has same sample order as DNEAobj object
T1Dmeta <- T1Dmeta[sampleNames(dnw), ]

#add new metadata to DNEAobj object
dnw <- includeMetadata(object = dnw, type = "sample", metadata = T1Dmeta)

}
\seealso{
\code{\link{featureNames}},\code{\link{sampleNames}},
}
\author{
Christopher Patsalis
}
