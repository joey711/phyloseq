\name{sampleMap<-}
\alias{assign-sampleMap}
\alias{sampleMap<-}
\title{Assign a new sampleMap to \code{x}}
\usage{
  sampleMap(x) <- value
}
\arguments{
  \item{x}{(Required). \code{\link{phyloseq-class}}}

  \item{value}{(Required). \code{\link{sampleMap-class}}}
}
\description{
  Assign a new sampleMap to \code{x}
}
\examples{
#
# data(ex1)
# # An example of pruning to just the first 10 samples in ex1
# ex2a <- prune_samples(sample.names(ex1)[1:10], ex1)
# # The following 3 lines produces an ex2b that is equal to ex2a
# ex2b <- ex1
# SD <- sampleMap(ex1)[1:10, ]
# sampleMap(ex2b) <- SD
# identical(ex2a, ex2b)
# print(ex2b)
}

