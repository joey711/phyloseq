\name{getslots.phyloseq}
\alias{getslots.phyloseq}
\title{Return the non-empty slot names of a phyloseq object.}
\usage{
  getslots.phyloseq(x)
}
\arguments{
  \item{x}{A \code{\link{phyloseq-class}} object. If
  \code{x} is a component data class, then just returns the
  class of \code{x}.}
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
  Like \code{\link{getSlots}}, but returns the class name
  if argument is component data object.
}
\examples{
#
}
\seealso{
  merge_phyloseq
}

