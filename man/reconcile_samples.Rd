\name{reconcile_samples}
\alias{reconcile_samples}
\title{Keep only sample-indices common to all components.}
\usage{
  reconcile_samples(x)
}
\arguments{
  \item{x}{An instance of phyloseq-class that contains 2 or
  more component data tables that in-turn describe
  samples.}
}
\value{
  A trimmed version of the argument, \code{x}, in which
  each component describes exactly the same set of samples.
  Class of \code{x} should be unchanged.
}
\description{
  This function is used internally as part of the
  infrastructure to ensure that component data types in a
  phyloseq-object describe exactly the same samples. In
  expected cases, a user will not need to invoke this
  function, because phyloseq objects are reconciled during
  instantiation by default.
}
\examples{
#
## data(ex1)
## reconcile_samples(ex1)
}
\seealso{
  \code{\link{reconcile_species}}
}
\keyword{internal}

