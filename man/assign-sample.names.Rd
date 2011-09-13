\name{sample.names<-}
\alias{assign-sample.names}
\alias{sample.names<-}
\title{Assign sample names to an object.}
\description{
  A convenience / syntax wrapper around
  \code{\link{prune_samples}}. This is typically used for
  subsetting an object to just those samples indicated by
  the character vector \code{value}.
}
\examples{
## data(ex1)
# # # ex3 <- ex1
# # # sample.names(ex3) <- sample(sample.names(ex1), 5)
# # # sample.names(ex3) <- sample.names(ex3)[sampleMap(ex3)[,"Gender"]=="A"]
}
\seealso{
  \code{\link{prune_samples}}
}

