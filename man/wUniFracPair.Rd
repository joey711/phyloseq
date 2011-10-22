\name{wUniFracPair}
\alias{wUniFracPair}
\title{Calculate weighted UniFrac for a pair of samples from an OTU table.}
\usage{
  wUniFracPair(OTU, tree, A, B, normalized=TRUE)
}
\arguments{
  \item{OTU}{otuTable in samples-by-species orientation}

  \item{tree}{object of class \code{phylo}}

  \item{A}{single character string matching the first
  sample ID in the pair}

  \item{B}{single character string matching the second
  sample ID in the pair}

  \item{normalized}{Logical. Should the output be
  normalized such that values range from 0 to 1 independent
  of branch length values? Default is \code{TRUE}.}
}
\value{
  numeric
}
\description{
  Somewhat an internal function for \code{\link{wUniFrac}},
  as usually interested in the weighted-UniFrac distances
  between many different samples. Returns a single numeric
  value, between 0 and 1 if normalized
}
\examples{
#
}
\seealso{
  wUniFrac
}

