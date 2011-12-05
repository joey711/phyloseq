\name{subset_samples}
\alias{subset_samples}
\title{Subset samples by sampleMap expression}
\usage{
  subset_samples(x, ...)
}
\arguments{
  \item{x}{A \code{sampleMap}, or a \code{phyloseq-class}
  object with a \code{sampleMap}. If the \code{sampleMap}
  slot is missing in \code{x}, then \code{x} will be
  returned as-is, and a warning will be printed to screen.}

  \item{...}{The subsetting expression that should be
  applied to the \code{sampleMap}. This is passed on to
  \code{\link{subset}}, see its documentation for more
  details.}
}
\value{
  A subsetted object with the same class as \code{x}.
}
\description{
  This is a convenience wrapper around the
  \code{\link{subset}} function. It is intended to allow
  subsetting complex experimental objects with one function
  call. In the case of \code{subset_samples}, the
  subsetting will be based on an expression related to the
  columns and values within the sampleMap slot of \code{x}.
}
\examples{
## data(ex1)
## ex2 <- subset_samples(ex1, Gender=="A")
## ex2 <- subset_samples(ex1, Diet==1)
}
\seealso{
  \code{\link{subset_species}}
}

