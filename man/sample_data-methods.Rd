\docType{methods}
\name{sample_data}
\alias{samData}
\alias{sam_data}
\alias{sampleData}
\alias{sample_data}
\alias{sample_data,ANY-method}
\alias{sample_data,data.frame-method}
\title{Build or access sample_data.}
\usage{
  sample_data(object, errorIfNULL=TRUE)

  sam_data(object, errorIfNULL = TRUE)

  sampleData(object, errorIfNULL = TRUE)

  samData(object, errorIfNULL = TRUE)
}
\arguments{
  \item{object}{(Required). A
  \code{\link{data.frame-class}}, or a
  \code{\link{phyloseq-class}} object.}

  \item{errorIfNULL}{(Optional). Logical. Should the
  accessor stop with an error if the slot is empty
  (\code{NULL})? Default \code{TRUE}.}
}
\value{
  A \code{\link{sample_data-class}} object representing the
  sample variates of an experiment.
}
\description{
  This is the suggested method for both constructing and
  accessing a table of sample-level variables
  (\code{\link{sample_data-class}}), which in the
  \code{\link{phyloseq-package}} is represented as a
  special extension of the \code{\link{data.frame-class}}.
  When the argument is a data.frame, sample_data() will
  create a sample_data-class object. In this case, the rows
  should be named to match the \code{sample.names} of the
  other objects to which it will ultimately be paired.
  Alternatively, if the first argument is an
  experiment-level (\code{\link{phyloseq-class}}) object,
  then the corresponding \code{sample_data} is returned.
  Like other accessors (see See Also, below), the default
  behavior of this method is to stop with an error if
  \code{object} is a \code{phyloseq-class} but does not
  contain a \code{sample_data}.
}
\details{
  Note that the \code{sam_data()} and \code{sampleMap()}
  functions are provided for convenience and backward
  compatibility, respectively, but should provide the exact
  same behavior as \code{sample_data()}.
}
\examples{
#
# data(GlobalPatterns)
# sample_data(GlobalPatterns)
## shorter (convenience) wrapper of sample_data()
# sam_data(GlobalPatterns)
}
\seealso{
  \code{\link{phy_tree}}, \code{\link{tax_table}},
  \code{\link{otu_table}} \code{\link{phyloseq}},
  \code{\link{merge_phyloseq}}
}

