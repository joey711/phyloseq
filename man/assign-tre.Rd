\docType{methods}
\name{tre<-}
\alias{assign-tre}
\alias{tre<-}
\alias{tre<-,phyloseq,phylo-method}
\title{Assign to tre an object/value.}
\usage{
  tre(x) <- value
}
\arguments{
  \item{x}{(Required). \code{\link{phyloseq-class}}}

  \item{value}{(Required). \code{\link{phylo-class}}}
}
\description{
  Assign to tre an object/value.
}
\examples{
#
# data(ex1)
# # An example of pruning to just the first 100 taxa in ex1.
# ex2a <- prune_species(species.names(ex1)[1:100], ex1)
# # The following 3 lines produces an ex2b that is equal to ex2a
# ex2b <- ex1
# tree <- prune_species(species.names(ex1)[1:100], tre(ex1))
# tre(ex2b) <- tree
# identical(ex2a, ex2b)
# print(ex2b)
}

