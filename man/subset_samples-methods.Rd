\docType{methods}
\name{subset_samples}
\alias{subset_samples}
\title{Subset samples by sampleData expression}
\usage{
  subset_samples(physeq, ...)
}
\arguments{
  \item{physeq}{A \code{\link{sampleData-class}}, or a
  \code{\link{phyloseq-class}} object with a
  \code{sampleData}. If the \code{sampleData} slot is
  missing in \code{physeq}, then \code{physeq} will be
  returned as-is, and a warning will be printed to screen.}

  \item{...}{The subsetting expression that should be
  applied to the \code{sampleData}. This is passed on to
  \code{\link{subset}}, see its documentation for more
  details.}
}
\value{
  A subsetted object with the same class as \code{physeq}.
}
\description{
  This is a convenience wrapper around the
  \code{\link{subset}} function. It is intended to allow
  subsetting complex experimental objects with one function
  call. The subsetting will be based on an expression
  related to the columns and values within the sampleData.
}
\examples{
# data(ex1)
 # ex2 <- subset_samples(ex1, Gender=="A")
 # ex2 <- subset_samples(ex1, Diet==1)
 ### Here is an example comparing subset_samples with prune_samples...
 # B_only_sample_names <- sample.names(sampleData(ex1)[(sampleData(ex1)[, "Gender"]=="B"),])
 # ex2 <- prune_samples(B_only_sample_names, ex1)
 # ex3 <- subset_samples(ex1, Gender=="B")
 # ## This should be TRUE.
 # identical(ex2, ex3)
}
\seealso{
  \code{\link{subset_species}}
}

