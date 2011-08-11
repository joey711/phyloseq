\name{threshrankfun}
\alias{threshrankfun}
\title{A closure version of the \code{threshrank} function.}
\usage{
  threshrankfun(thresh, keep0s = FALSE, ...)
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
  A single-argument function with the options to threshrank
  set.
}
\description{
  Takes the same arguments as \code{\link{threshrank}},
  except for \code{x}, because the output is a function
  rather than a rank-transformed numeric. This approach is
  useful for creating an input to a higher-order function,
  like "filterfun", that require a single-argument function
  as input.
}
\examples{
#
}
\seealso{
  topk
}

