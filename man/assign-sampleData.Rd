\name{sampleData<-}
\alias{assign-sampleData}
\alias{samData<-}
\alias{sampleData<-}
\title{Assign a new sampleData to \code{x}}
\usage{
  sampleData(x) <- value

  samData(x) <- value
}
\arguments{
  \item{x}{(Required). \code{\link{phyloseq-class}}}

  \item{value}{(Required). \code{\link{sampleData-class}}}
}
\description{
  Assign a new sampleData to \code{x}
}
\examples{
#
# data(ex1)
# # An example of pruning to just the first 10 samples in ex1
# ex2a <- prune_samples(sample.names(ex1)[1:10], ex1)
# # The following 3 lines produces an ex2b that is equal to ex2a
# ex2b <- ex1
# SD <- sampleData(ex1)[1:10, ]
# sampleData(ex2b) <- SD
# identical(ex2a, ex2b)
# print(ex2b)
}

