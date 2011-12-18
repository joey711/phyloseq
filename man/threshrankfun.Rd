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
  \code{\link{transformsamplecounts}}.
}
\examples{
data(ex1)
## These three approaches result in identical otuTable
(x1 <- transformsamplecounts( otuTable(ex1), threshrankfun(500)) )
(x2 <- otuTable(apply(otuTable(ex1), 2, threshrankfun(500)), speciesAreRows(ex1)) )
identical(x1, x2)
(x3 <- otuTable(apply(otuTable(ex1), 2, threshrank, thresh=500), speciesAreRows(ex1)) )
identical(x1, x3)
}
\seealso{
  \code{\link{transformsamplecounts}},
  \code{\link{threshrankfun}}, \code{\link{threshrank}}
}

