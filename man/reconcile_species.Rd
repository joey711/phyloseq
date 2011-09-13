\name{reconcile_species}
\alias{reconcile_species}
\title{Keep only species-indices common to all components.}
\usage{
  reconcile_species(x)
}
\arguments{
  \item{x}{An object of the phyloseq package that contains
  2 or more components data tables that in-turn describe
  species/taxa. E.g. An otuTree object, or an otuTax
  object.}
}
\value{
  A trimmed version of the argument, \code{x}, in which
  each component describes exactly the same set of
  species/taxa. Class of \code{x} should be unchanged.
}
\description{
  This function is used internally as part of the
  infrastructure to ensure that component data types in a
  phyloseq-object have exactly the same taxa/species. It
  relies heavily on the \code{\link{prune_species}} S4
  methods to perform the actual trimming. In expected
  cases, a user will not need to invoke this function,
  because phyloseq objects are reconciled during
  instantiation by default.
}
\examples{
#
## data(ex1)
## head(intersect_species(ex1), 10)
## reconcile_species(ex1)
}
\seealso{
  \code{\link{reconcile_samples}}, \code{\link{Reduce}}
}

