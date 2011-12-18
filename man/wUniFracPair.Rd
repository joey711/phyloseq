\name{wUniFracPair}
\alias{wUniFracPair}
\title{Calculate weighted UniFrac for a pair of samples from an abundance matrix.}
\usage{
  wUniFracPair(OTU, tree, A, B, normalized=TRUE)
}
\arguments{
  \item{OTU}{(Required). An abundance matrix in
  samples-by-species orientation.}

  \item{tree}{(Required). Object of class
  \code{\link[ape]{phylo}}}

  \item{A}{(Required). single character string matching the
  first sample ID in the pair}

  \item{B}{(Required). single character string matching the
  second sample ID in the pair}

  \item{normalized}{(Optional). Logical. Should the output
  be normalized such that values range from 0 to 1
  independent of branch length values? Default is
  \code{TRUE}.}
}
\value{
  A single number between 0, 1.
}
\description{
  Somewhat an internal function for \code{\link{UniFrac}},
  as usually interested in the weighted-UniFrac distances
  between many different samples. Returns a single numeric
  value, between 0 and 1 if normalized
}
\seealso{
  \code{\link{UniFrac}}
}
\keyword{internal}

