\name{splat.phyloseq.objects}
\alias{splat.phyloseq.objects}
\title{Convert phyloseq-class into a named list of its non-empty components.}
\usage{
  splat.phyloseq.objects(x)
}
\arguments{
  \item{x}{An object of a class defined by the
  phyloseq-package. Component data and complex classes are
  both acceptable.}
}
\value{
  A named list, where each element is a component object
  that was contained in the argument, \code{x}. Each
  element is named for the object class it contains If
  \code{x} is already a component data object, then a list
  of length (1) is returned, also named.
}
\description{
  Convert phyloseq-class into a named list of its non-empty
  components.
}
\examples{
#
}
\seealso{
  merge_phyloseq
}

