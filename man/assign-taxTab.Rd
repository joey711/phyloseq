\name{taxTab<-}
\alias{assign-taxTab}
\alias{taxTab<-}
\title{Assign a new Taxonomy Table to \code{x}}
\usage{
  taxTab(x) <- value
}
\arguments{
  \item{x}{(Required). \code{\link{phyloseq-class}}}

  \item{value}{(Required).
  \code{\link{taxonomyTable-class}}}
}
\description{
  Assign a new Taxonomy Table to \code{x}
}
\examples{
#
# data(ex1)
# # An example of pruning to just the first 100 taxa in ex1.
# ex2a <- prune_species(species.names(ex1)[1:100], ex1)
# # The following 3 lines produces an ex2b that is equal to ex2a
# ex2b <- ex1
# TT <- taxTab(ex1)[1:100, ]
# taxTab(ex2b) <- TT
# identical(ex2a, ex2b)
# print(ex2b)
}

