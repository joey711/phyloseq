\name{intersect_species}
\alias{intersect_species}
\title{Returns the intersection of species for the components of x}
\usage{
  intersect_species(x)
}
\arguments{
  \item{x}{(Required). A \code{\link{phyloseq-class}}
  object that contains 2 or more components that in-turn
  describe species/taxa.}
}
\value{
  Returns a character vector of only those species that are
  present in all species-describing components of \code{x}.
}
\description{
  This function is used internally as part of the
  infrastructure to ensure that component data types in a
  phyloseq-object have exactly the same taxa/species. It
  relies heavily on the \code{\link{Reduce}} function to
  determine the strictly common species.
}
\examples{
#
## data(GlobalPatterns)
## head(intersect_species(GlobalPatterns), 10)
}
\seealso{
  \code{\link{reconcile_species}}, \code{\link{Reduce}}
}
\keyword{internal}

