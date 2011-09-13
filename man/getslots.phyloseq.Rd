\name{getslots.phyloseq}
\alias{getslots.phyloseq}
\title{return the slot names of phyloseq objects.}
\usage{
  getslots.phyloseq(x)
}
\arguments{
  \item{x}{An object of a class defined by the
  phyloseq-package. Component data and complex classes are
  both acceptable.}
}
\value{
  identical to getSlots. A named character vector of the
  slot classes of a particular S4 class, where each element
  is named by the slot name it represents. If \code{x} is a
  component data object, then a vector of length (1) is
  returned, named according to its slot name in the
  higher-order objects.
}
\description{
  Like getSlots, but returns the class name if argument is
  component data object.
}
\examples{
#
}
\seealso{
  merge_phyloseq
}

