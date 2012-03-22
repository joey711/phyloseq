\docType{methods}
\name{tre<-}
\alias{assign-tre}
\alias{tre<-}
\alias{tre<-,phyloseq,phylo-method}
\alias{tre<-,phyloseq,phyloseq-method}
\title{Assign a (new) phylogenetic tree to \code{x}}
\usage{
  tre(x) <- value
}
\arguments{
  \item{x}{(Required). \code{\link{phyloseq-class}}}

  \item{value}{(Required). \code{\link{phylo-class}}, or
  \code{\link{phyloseq-class}}}
}
\description{
  Assign a (new) phylogenetic tree to \code{x}
}
\examples{
#
# data(GlobalPatterns)
# # An example of pruning to just the first 100 taxa in GlobalPatterns.
# ex2a <- prune_species(species.names(GlobalPatterns)[1:100], GlobalPatterns)
# # The following 3 lines produces an ex2b that is equal to ex2a
# ex2b <- GlobalPatterns
# tree <- prune_species(species.names(GlobalPatterns)[1:100], tre(GlobalPatterns))
# tre(ex2b) <- tree
# identical(ex2a, ex2b)
# print(ex2b)
# # Example adding a phylo tree from phyloseq class
# ex2c <- phyloseq(otuTable(ex2b), sampleData(ex2b), taxTab(ex2b))
# tre(ex2c) <- ex2b
# identical(ex2b, ex2c)
}

