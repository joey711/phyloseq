\docType{methods}
\name{otuTable}
\alias{otuTable}
\alias{otuTable,data.frame-method}
\alias{otuTable,matrix-method}
\alias{otuTable,otuTable-method}
\alias{otuTable,phyloseq-method}
\title{Build or access the otuTable.}
\usage{
  otuTable(object, speciesAreRows, errorIfNULL=TRUE)
}
\arguments{
  \item{object}{(Required). An integer matrix,
  \code{\link{otuTable-class}}, or
  \code{\link{phyloseq-class}}.}

  \item{speciesAreRows}{(Conditionally optional). Logical;
  of length 1. Ignored unless \code{object} is a matrix, in
  which case it is is required.}

  \item{errorIfNULL}{(Optional). Logical. Should the
  accessor stop with an error if the slot is empty
  (\code{NULL})? Default \code{TRUE}. Ignored if
  \code{object} argument is a matrix (constructor invoked
  instead).}
}
\value{
  An \code{\link{otuTable-class}} object.
}
\description{
  This is the suggested method for both constructing and
  accessing Operational Taxonomic Unit (OTU) abundance
  (\code{\link{otuTable-class}}) objects. When the first
  argument is a matrix, otuTable() will attempt to create
  and return an otuTable-class object, which further
  depends on whether or not \code{speciesAreRows} is
  provided as an additional argument. Alternatively, if the
  first argument is an experiment-level
  (\code{\link{phyloseq-class}}) object, then the
  corresponding \code{otuTable} is returned.
}
\examples{
#
# data(GlobalPatterns)
# otuTable(GlobalPatterns)
}
\seealso{
  \code{\link{tre}}, \code{\link{sampleData}},
  \code{\link{taxTab}} \code{\link{phyloseq}},
  \code{\link{merge_phyloseq}}
}

