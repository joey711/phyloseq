\docType{methods}
\name{prune_samples}
\alias{prune_samples}
\alias{prune_samples,character,otuTable-method}
\alias{prune_samples,character,phyloseq-method}
\alias{prune_samples,character,sampleData-method}
\title{Prune unwanted samples from a phyloseq object.}
\usage{
  prune_samples(samples, x)
}
\arguments{
  \item{samples}{A character vector of the samples in
  object x that you want to keep.}

  \item{x}{A phyloseq object.}
}
\value{
  The class of the object returned by \code{prune_samples}
  matches the class of the phyloseq object, \code{x}.
}
\description{
  An S4 Generic method for removing (pruning) unwanted
  samples.
}
\examples{
#
# data(GlobalPatterns)
# GP <- GlobalPatterns
# B_only_sample_names <- sample.names(sampleData(GP)[(sampleData(GP)[, "Gender"]=="B"),])
# ex2 <- prune_samples(B_only_sample_names, GP)
# ex3 <- subset_samples(GP, Gender=="B")
# ## This should be TRUE.
# identical(ex2, ex3)
# ## Here is a simpler example: Make new object with only the first 5 samples
# ex4 <- prune_samples(sample.names(GP)[1:5], GP)
}
\seealso{
  \code{\link{subset_samples}}
}

