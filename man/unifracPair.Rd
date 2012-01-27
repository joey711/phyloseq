\name{unifracPair}
\alias{unifracPair}
\title{Calculate (unweighted) UniFrac for a pair of samples from an occurrence matrix.}
\usage{
  unifracPair(occ, tree, A, B)
}
\arguments{
  \item{occ}{(Required). A matrix. A samples-by-species
  occurrence matrix.}

  \item{tree}{(Required). A rooted phylogenetic tree.
  \code{phylo} class.}

  \item{A}{(Required). Single character string. The name of
  sample \code{"A"}.}

  \item{B}{(Required). Single character string. The name of
  sample \code{"B"}.}
}
\value{
  A single number between 0, 1.
}
\description{
  Somewhat an internal function for \code{\link{UniFrac}},
  as an investigator is usually interested in the UniFrac
  distances between many different samples. Returns a
  single numeric value, between 0 and 1.
}
\seealso{
  See the main function, \code{\link{UniFrac}}.
}
\keyword{internal}

