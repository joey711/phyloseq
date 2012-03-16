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
# data(GlobalPatterns)
# # An example of pruning to just the first 100 taxa in GlobalPatterns.
# ex2a <- prune_species(species.names(GlobalPatterns)[1:100], GlobalPatterns)
# # The following 3 lines produces an ex2b that is equal to ex2a
# ex2b <- GlobalPatterns
# TT <- taxTab(GlobalPatterns)[1:100, ]
# taxTab(ex2b) <- TT
# identical(ex2a, ex2b)
# print(ex2b)
}

