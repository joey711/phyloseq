\name{threshrankfun}
\alias{threshrankfun}
\title{A closure version of the \code{threshrank} function.}
\usage{
  threshrankfun(thresh, keep0s=FALSE, ...)
}
\arguments{
  \item{thresh}{A single numeric value giving the
  threshold.}

  \item{keep0s}{A logical determining whether 0's in
  \code{x} should remain a zero-value in the output. If
  FALSE, zeros are treated as any other value.}

  \item{...}{Further arguments passes to the
  \code{\link{rank}} function.}
}
\value{
  A single-argument function with the options to
  \code{\link{threshrank}} set.
}
\description{
  Takes the same arguments as \code{\link{threshrank}},
  except for \code{x}, because the output is a
  single-argument function rather than a rank-transformed
  numeric. This is useful for higher-order functions that
  require a single-argument function as input, like
  \code{\link{transform_sample_counts}}.
}
\examples{
data(GlobalPatterns)
GP <- GlobalPatterns
## These three approaches result in identical otu_table
(x1 <- transform_sample_counts( otu_table(GP), threshrankfun(500)) )
(x2 <- otu_table(apply(otu_table(GP), 2, threshrankfun(500)), taxa_are_rows(GP)) )
identical(x1, x2)
(x3 <- otu_table(apply(otu_table(GP), 2, threshrank, thresh=500), taxa_are_rows(GP)) )
identical(x1, x3)
}
\seealso{
  \code{\link{transform_sample_counts}},
  \code{\link{threshrankfun}}, \code{\link{threshrank}}
}

