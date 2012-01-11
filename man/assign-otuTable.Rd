\docType{methods}
\name{otuTable<-}
\alias{assign-otuTable}
\alias{otuTable<-}
\alias{otuTable<-,otuTable,otuTable-method}
\alias{otuTable<-,phyloseq,otuTable-method}
\title{Assign a new OTU Table to \code{x}}
\arguments{
  \item{x}{(Required). \code{\link{phyloseq-class}}}

  \item{value}{(Required). \code{\link{otuTable-class}}}
}
\description{
  Assign a new OTU Table to \code{x}
}
\examples{
# data(ex1)
# # An example of pruning to just the first 100 taxa in ex1.
# ex2a <- prune_species(species.names(ex1)[1:100], ex1)
# # The following 3 lines produces an ex2b that is equal to ex2a
# ex2b <- ex1
# OTU <- otuTable(ex1)[1:100, ]
# otuTable(ex2b) <- OTU
# identical(ex2a, ex2b)
# print(ex2b)
}

