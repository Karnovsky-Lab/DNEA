% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities-exported.R
\name{addExpressionData}
\alias{addExpressionData}
\title{Include custom normalized data in the DNEAobj object}
\usage{
addExpressionData(object, data)
}
\arguments{
\item{object}{A \code{\link{DNEAobj}} object}

\item{data}{An \emph{m x n} numeric matrix of custom-normalized
expression expression data}
}
\value{
A \code{\link{DNEAobj}} object with the added expression data in
the @assays slot
}
\description{
This function allows the user to input custom-normalized data into the
DNEAobj object for use in DNEA analysis.
}
\examples{
#import example data
data(TEDDY)
data(dnw)
data(T1Dmeta)

#transpose TEDDY data
TEDDY <- t(TEDDY)

#make sure metadata and expression data are in same order
T1Dmeta <- T1Dmeta[rownames(TEDDY),]

dat <- list('DM:control' = TEDDY[T1Dmeta$group == "DM:control",],
            'DM:case' = TEDDY[T1Dmeta$group == "DM:case",])

#log-transform and median center the expression data without scaling
newdat <- list()
for(cond in seq(length(dat))){

  group_dat <- dat[[cond]]
  for(i in seq(1, ncol(group_dat))){
    metab_median = median(group_dat[, i], na.rm = TRUE)
    metab_range = range(group_dat[, i], na.rm = TRUE)
    scale_factor = max(abs(metab_range - metab_median))
    group_dat[, i] <- (group_dat[, i] - metab_median) / scale_factor

    rm(metab_median, metab_range, scale_factor)
  }

  group_dat <- group_dat[rownames(dat[[cond]]),colnames(dat[[cond]])]
  group_dat <- t(group_dat)
  newdat <- append(newdat, list(group_dat))

  rm(i, group_dat)
}

#add names
names(newdat) <- names(dat)

#add data
dnw <- addExpressionData(object = dnw, data = newdat)

}
\seealso{
\code{\link{createDNEAobject}},\code{\link{DNEAobj}},
}
\author{
Christopher Patsalis
}
