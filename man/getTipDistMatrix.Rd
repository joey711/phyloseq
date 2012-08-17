\name{getTipDistMatrix}
\alias{getTipDistMatrix}
\alias{gettipdistmatrix}
\title{An internal wrapper function on \code{\link[ape]{cophenetic.phylo}}}
\arguments{
  \item{tree}{\code{phylo}}

  \item{byRootFraction}{Should the distance be calculated
  according to fractional distance to the root? If
  \code{FALSE}, the distance is instead the patristic
  distance as calculated by cophenetic.phylo. Default
  \code{FALSE}.}
}
\value{
  character matrix. First column is the complete match,
  followed by one for each capture group
}
\description{
  This is useful for determining tips that should be
  combined.
}
\seealso{
  tip_glom
}
\keyword{internal}

