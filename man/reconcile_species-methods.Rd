\name{reconcile_species}
\alias{reconcile_species}
\alias{reconcile_species,phyloseq-method}
\title{Keep only species-indices common to all components.}
\usage{
  reconcile_species(x)
}
\arguments{
  \item{x}{(Required). A phyloseq object. Only meaningful
  if \code{x} has at least two non-empty slots of the
  following slots that describe species:
  \code{\link{otuTable}}, \code{\link{taxTab}},
  \code{\link{tre}}.}
}
\value{
  A trimmed version of the argument, \code{x}, in which
  each component describes exactly the same set of
  species/taxa. Class of return should match argument,
  \code{x}.
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
## data(GlobalPatterns)
## head(intersect_species(GlobalPatterns), 10)
## reconcile_species(GlobalPatterns)
# # data(phylocom)
# # tree <- phylocom$phylo
# # OTU  <- otuTable(phylocom$sample, speciesAreRows=FALSE)
# # ex3  <- phyloseq(OTU, tree)
# # reconcile_species(ex3)
# # intersect_species(ex3)
}
\seealso{
  \code{\link{reconcile_samples}}, \code{\link{Reduce}}
}
\keyword{internal}

