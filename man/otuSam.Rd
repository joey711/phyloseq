\name{otuSam}
\alias{otusam}
\alias{otuSam}
\title{Create a subset object with just the \code{otuTable} and \code{sampleMap} slots.}
\usage{
  otuSam(object)
}
\arguments{
  \item{object}{An object among the set of classes defined
  by the phyloseq package that contain both a sampleMap and
  an otuTable.}
}
\value{
  A phyloseq-package object of class \code{otuSam}.
}
\description{
  Create a subset object with just the \code{otuTable} and
  \code{sampleMap} slots.
}
\examples{
#
## data(ex1)
## class(ex1)
## otuSam(ex1)
}

