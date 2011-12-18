\name{threshrank}
\alias{threshrank}
\title{Thresholded rank transformation.}
\usage{
  threshrank(x, thresh, keep0s=FALSE, ...)
}
\arguments{
  \item{x}{The numeric vector to transform}

  \item{thresh}{A single numeric value giving the
  threshold.}

  \item{keep0s}{A logical determining whether 0's in
  \code{x} should remain a zero-value in the output. If
  FALSE, zeros are treated as any other value.}

  \item{...}{Further arguments passes to the
  \code{\link{rank}} function.}
}
\value{
  A ranked, (optionally) thresholded numeric vector with
  length equal to \code{x}. Default arguments to
  \code{rank} are used, unless provided as additional
  arguments.
}
\description{
  The lowest \code{thresh} values in \code{x} all get the
  value 'thresh'.
}
\examples{
#
threshrank(sample(0:10, 100, TRUE), 50, keep0s=TRUE)
}
\seealso{
  \code{\link{rank}}, \code{\link{threshrankfun}}
}

