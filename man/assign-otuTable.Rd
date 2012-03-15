\docType{methods}
\name{otuTable<-}
\alias{assign-otuTable}
\alias{otuTable<-}
\alias{otuTable<-,otuTable,otuTable-method}
\alias{otuTable<-,phyloseq,otuTable-method}
\title{Assign a new OTU Table to \code{x}}
\usage{
  otuTable(x) <- value
}
\arguments{
  \item{x}{(Required). \code{\link{phyloseq-class}}}

  \item{value}{(Required). \code{\link{otuTable-class}}}
}
\description{
  Assign a new OTU Table to \code{x}
}
\examples{
# data(GlobalPatterns)
# # An example of pruning to just the first 100 taxa in GlobalPatterns.
# ex2a <- prune_species(species.names(GlobalPatterns)[1:100], GlobalPatterns)
# # The following 3 lines produces an ex2b that is equal to ex2a
# ex2b <- GlobalPatterns
# OTU <- otuTable(GlobalPatterns)[1:100, ]
# otuTable(ex2b) <- OTU
# identical(ex2a, ex2b)
# print(ex2b)
}

