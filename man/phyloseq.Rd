\name{phyloseq}
\alias{phyloseq}
\title{Build objects for the phyloseq package.}
\usage{
  phyloseq(...)
}
\arguments{
  \item{...}{One or more component objects among the set of
  classes defined by the phyloseq package, as well as
  \code{phylo4} objects (defined by the phylobase package).
  Each argument should be a different class.  For combining
  multiple components of the same class, or multiple
  higher-order classes, use the
  \code{\link{merge_phyloseq}} function. Unlike in earlier
  versions, the arguments to phyloseq do not need to be
  named, and the order of the arguments does not matter.}
}
\value{
  The class of the returned object depends on the argument
  class(es). To construct a H.O. object, two or more
  component data objects must be provided in the argument
  list.  Otherwise, the order of arguments does not matter.
  If a single component-class object is provided, it is
  simply returned as-is.
}
\description{
  \code{phyloseq()} is a constructor method, This is the
  main method suggested for constructing \code{phyloseq}
  higher-order objects from their components.
}
\examples{
#
# # data(ex1)
# # phyloseq(sampleMap(ex1), otuTable(ex1))
# # phyloseq(otuTable(ex1), tre(ex1))
# # phyloseq(taxTab(ex1), otuTable(ex1))
# # phyloseq(tre(ex1), otuTable(ex1), sampleMap(ex1))
# # phyloseq(otuTable(ex1), taxTab(ex1), sampleMap(ex1))
# # phyloseq(otuTable(ex1), tre(ex1), taxTab(ex1), sampleMap(ex1))
}
\seealso{
  \code{\link{merge_phyloseq}}
}

