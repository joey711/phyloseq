\name{topp}
\alias{topp}
\title{The most abundant \code{p} fraction of taxa}
\usage{
  topp(p, na.rm = TRUE)
}
\arguments{
  \item{p}{A numeric of length 1, indicating what fraction
  of the most abundant taxa should be kept.}

  \item{na.rm}{A logical. Should \code{NA}s be removed.
  Default is \code{TRUE}.}
}
\value{
  A function (enclosure), suitable for
  \code{\link{filterfunSample}}, that will return
  \code{TRUE} for each element in the most abundant p
  fraction of taxa.
}
\description{
  The most abundant \code{p} fraction of taxa
}

