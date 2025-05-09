% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initiator.R
\name{massDataset2DNEA}
\alias{massDataset2DNEA}
\title{Initialize a DNEA object from a mass_dataset object}
\usage{
massDataset2DNEA(project_name, object, group_label_col, scaled_input = FALSE)
}
\arguments{
\item{project_name}{A character string name for the experiment.}

\item{object}{a mass_dataset object.}

\item{group_label_col}{A character string corresponding to the
column in the sample metadata stored in the mass_dataset object
to use as the group labels.}

\item{scaled_input}{A TRUE/FALSE indicating whether the input data
is already normalized}
}
\value{
A \code{\link[=DNEA-class]{DNEA}} object.
}
\description{
This function takes as input a
mass_dataset-class object from the massdataset package to
initiate a \code{\link[=DNEA-class]{DNEA}} object. Differential
expression analysis is performed using a student's T-test
and Benjamini-Hochberg for multiple-testing corrections.
Diagnostic testing is done on the input data by checking
the minimum eigen value and condition number of the
expression data for each experimental condition.
\emph{\strong{NOTE: the massdataset package from the
tidymass software suite must be installed to use
this function. Please see
\url{https://massdataset.tidymass.org/} for more installation
instructions}}
\subsection{IMPORTANT:}{

Special attention should be given to the diagnostic criteria that is
output. The minimum eigen value and condition number are calculated for
the whole data set as well as for each condition to determine mathematic
stability of the data set and subsequent results from a GGM model. More
information about interpretation can be found in the
\strong{\emph{Details}} section below.
}
}
\details{
\subsection{Diagnostics Motivation}{

Negative or zero eigenvalues in a data set can represent
instability in that portion of the matrix, thereby invalidating
parametric statistical methods and creating unreliable results. In this
function, the minimum eigenvalue of the data set is calculated by first
creating a pearson correlation matrix of the data. Instability may then
occur for a number of reasons, but one common cause is highly correlated
features (in the positive and negative direction). \cr

Regularization often takes care of this problem by arbitrarily
selecting one of the variables in a highly correlated group and removing
the rest. We have developed DNEA to be very robust in situations where
\strong{\emph{p >> n}} by optimizing the model via several regularization
steps (\emph{please see} \code{\link{BICtune}} \emph{and}
\code{\link{stabilitySelection}}) that may handle such problems without
intervention, however, the user can also pre-emptively collapse
highly-correlated features into a single group via
\code{\link{aggregateFeatures}}.
}

\subsection{Benefits of Feature Aggregation}{

When your dataset contains highly correlated features, we recommend
aggregating features into related groups - such as highly-correlated
features of a given class of molecules (ie. many fatty acids,
carnitines, etc.) - because the user then has more control over which
variables are included in the model. Without collapsing, the model
regularization may result in one of the features within a class being
included and some or all of the remaining features being removed. By
collapsing first, you retain the signal from all of the features in the
collapsed group and also have information pertaining to which features
are highly correlated and will therefore have similar
feature-feature associations.
}
}
\examples{
#load data
data(TEDDY)
data(T1Dmeta)
data(metab_data)

#make sure metadata and expression data are in same order
T1Dmeta <- T1Dmeta[colnames(TEDDY),]
T1Dmeta <- T1Dmeta[, c(6,7,7)]
colnames(T1Dmeta) <- c("sample_id", "group", "class")

metab_data <- metab_data[rownames(TEDDY), ]

sample_info_note = data.frame(name = c("sample_id", "group", "class"),
                              meaning = c("sample", "group", "class"))
variable_info_note = data.frame(name = c("variable_id", "mz", "rt"),
                                meaning = c("variable_id", "mz", "rt"))
if (require(massdataset)) {
#create mass_dataset object from TEDDY
object <- massdataset::create_mass_dataset(expression_data = data.frame(TEDDY),
                                           sample_info = T1Dmeta,
                                           variable_info = metab_data,
                                           sample_info_note = sample_info_note,
                                           variable_info_note = variable_info_note)

DNEA <- massDataset2DNEA(project_name = "mass_dataset",
                             object = object,
                             group_label_col = "group")
}

}
\seealso{
\code{\link{BICtune}}, \code{\link{stabilitySelection}},
\code{\link{createDNEAobject}}
}
\author{
Christopher Patsalis
}
