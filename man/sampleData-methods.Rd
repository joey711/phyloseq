\docType{methods}
\name{sampleData}
\alias{samData}
\alias{sampleData}
\alias{sampleData,ANY-method}
\alias{sampleData,data.frame-method}
\alias{sampleMap}
\title{Build or access sampleData.}
\usage{
  samData(object, errorIfNULL=TRUE)

  sampleMap(object, errorIfNULL=TRUE)
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
  A \code{\link{sampleData-class}} object representing the
  sample variates of an experiment.
}
\description{
  This is the suggested method for both constructing and
  accessing a table of sample-level variables
  (\code{\link{sampleData-class}}), which in the
  \code{\link{phyloseq-package}} is represented as a
  special extension of the \code{\link{data.frame-class}}.
  When the argument is a data.frame, sampleData() will
  create a sampleData-class object. In this case, the rows
  should be named to match the \code{sample.names} of the
  other objects to which it will ultimately be paired.
  Alternatively, if the first argument is an
  experiment-level (\code{\link{phyloseq-class}}) object,
  then the corresponding \code{sampleData} is returned.
  Like other accessors (see See Also, below), the default
  behavior of this method is to stop with an error if
  \code{object} is a \code{phyloseq-class} but does not
  contain a \code{sampleData}.
}
\details{
  Note that the \code{samData()} and \code{sampleMap()}
  functions are provided for convenience and backward
  compatibility, respectively, but should provide the exact
  same behavior as \code{sampleData()}.
}
\examples{
#
# data(ex1)
# sampleData(ex1)
## shorter (convenience) wrapper of sampleData()
# samData(ex1)
}
\seealso{
  \code{\link{tre}}, \code{\link{taxTab}},
  \code{\link{otuTable}} \code{\link{phyloseq}},
  \code{\link{merge_phyloseq}}
}

