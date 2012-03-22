\name{taxTab<-}
\alias{assign-taxTab}
\alias{taxTab<-}
\title{Assign a (new) Taxonomy Table to \code{x}}
\usage{
  taxTab(x) <- value
}
\arguments{
  \item{x}{(Required). \code{\link{phyloseq-class}}}

  \item{value}{(Required).
  \code{\link{taxonomyTable-class}}. Alternatively,
  \code{value} can be a \code{\link{phyloseq-class}} that
  has a \code{\link{taxTab}} component, or a
  \code{\link{matrix-class}} that can be coerced to a
  \code{\link{taxonomyTable-class}} with row indices that
  match at least some of the \code{\link{species.names}} of
  \code{x}.}
}
\description{
  Assign a (new) Taxonomy Table to \code{x}
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
# # 2 examples adding a taxTab component from phyloseq or matrix classes
# ex2c <- phyloseq(otuTable(ex2b), sampleData(ex2b), tre(ex2b))
# taxTab(ex2c) <- ex2b
# identical(ex2a, ex2c)
# ex2c <- phyloseq(otuTable(ex2b), sampleData(ex2b), tre(ex2b))
# taxTab(ex2c) <- as(taxTab(ex2b), "matrix")
# identical(ex2a, ex2c)
}

