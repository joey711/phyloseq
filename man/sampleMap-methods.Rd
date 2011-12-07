\docType{methods}
\name{sampleMap}
\alias{samplemap}
\alias{sampleMap}
\alias{sampleMap,ANY-method}
\alias{sampleMap,data.frame-method}
\title{Build or access the sampleMap.}
\usage{
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
  A \code{\link{sampleMap-class}} object representing the
  sample variates of an experiment.
}
\description{
  This is the suggested method for both constructing and
  accessing a table of sample-level variables
  (\code{\link{sampleMap-class}}), which in the
  \code{\link{phyloseq-package}} is represented as a
  special extension of the \code{\link{data.frame-class}}.
  When the argument is a data.frame, sampleMap() will
  create a sampleMap-class object. In this case, the rows
  should be named to match the \code{sample.names} of the
  other objects to which it will ultimately be paired.
  Alternatively, if the first argument is an
  experiment-level (\code{\link{phyloseq-class}}) object,
  then the corresponding \code{sampleMap} is returned. Like
  other accessors (see See Also, below), the default
  behavior of this method is to stop with an error if
  \code{object} is a \code{phyloseq-class} but does not
  contain a \code{sampleMap}.
}
\examples{
#
# data(ex1)
# sampleMap(ex1)
}
\seealso{
  \code{\link{tre}}, \code{\link{taxTab}},
  \code{\link{otuTable}} \code{\link{phyloseq}},
  \code{\link{merge_phyloseq}}
}

