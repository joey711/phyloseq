\name{UFwi}
\alias{UFwi}
\title{Internal function for \code{wUniFrac}.}
\usage{
  UFwi(edge, samples, OTU, tree, AT = sum(OTU[samples[1],
  ]), BT = sum(OTU[samples[2], ]))
}
\arguments{
  \item{edge}{The edge index}

  \item{samples}{Character vector of length 2, giving the
  pair of samples under comparison.}

  \item{OTU}{\code{otuTable} object in samples-by-species
  orientation}

  \item{tree}{object of class \code{phylo}}
}
\value{
  character vector
}
\description{
  A function that takes a phylo object (tree) and an
  edge-index as input, and returns the edge-weight term for
  weighted UniFrac.
}
\examples{
#
}
\seealso{
  wUniFrac
}
\keyword{internal}
\keyword{internal2tips}

