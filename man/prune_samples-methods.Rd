\docType{methods}
\name{prune_samples}
\alias{prune_samples}
\alias{prune_samples,character,otuTable-method}
\alias{prune_samples,character,phyloseq-method}
\alias{prune_samples,character,sampleMap-method}
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
}

