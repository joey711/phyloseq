\name{otuTax}
\alias{otutax}
\alias{otuTax}
\title{Create a subset object with just the \code{otuTable} and \code{taxTab} slots.}
\usage{
  otuTax(object)
}
\arguments{
  \item{object}{An object among the set of classes defined
  by the phyloseq package that contain a taxonomyTable, and
  an otuTable.}
}
\value{
  A phyloseq-package object of class \code{otuSamTax}.
}
\description{
  This is a convenience wrapper of \code{\link{phyloseq}}
  that constructs an otuTax-class object from a (equal or)
  more complex object that contains an otuTable and taxTab
  slots, effectively creating a subsetted object.
}
\details{
  This is an accessor method only. For building an otuTax
  object from its components, the \code{\link{phyloseq}}
  function is suggested.
}
\examples{
#
}

