\name{species.names<-}
\alias{assign-species.names}
\alias{species.names<-}
\title{Assign species names to an object.}
\description{
  A convenience / syntax wrapper around
  \code{\link{prune_samples}}. This is typically used for
  subsetting an object to just those species indicated by
  the character vector \code{value}.
}
\seealso{
  prune_species # # # data(ex1) ### subset ex1 to just 5
  random species # # # species.names(ex1) <-
  sample(species.names(ex1), 5) ### remove those species
  with only 1 or fewer individuals in ex1. # # #
  species.names(ex1) <- species.names(ex1)[speciesSums(ex1)
  > 2]
}

