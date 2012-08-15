\name{ufnum}
\alias{ufnum}
\title{Internal function for unweighted UniFrac edge-weight.}
\usage{
  ufnum(edge, samples, occ, tree)
}
\arguments{
  \item{edge}{The edge index}

  \item{samples}{Character vector of length 2, giving the
  pair of samples under comparison.}

  \item{occ}{\code{otu_table} object in samples-by-species
  orientation}

  \item{tree}{object of class \code{phylo}}
}
\value{
  character vector
}
\description{
  A function that takes a phylo object (tree) and an
  edge-index as input, and returns the edge-weight term for
  unweighted UniFrac.
}
\seealso{
  UniFrac
}
\keyword{internal}

