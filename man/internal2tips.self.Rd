\name{internal2tips.self}
\alias{internal2tips.self}
\title{Internal function for \code{\link{wUniFrac}}.}
\usage{
  internal2tips.self(phy, int.node, return.names = FALSE)
}
\arguments{
  \item{phy}{object of class \code{phylo}.}

  \item{int.node}{the internal node you want to get the
  descendants of.}
}
\value{
  character vector
}
\description{
  A modified version of the \code{\link{internal2tips}}
  function, such that when a tip is provided as int.node,
  that tip is returned. This is a more intuitive behavior
  than the picante version, which returns NULL, and
  currently used in the \code{\link{wUniFrac}}.
}
\examples{
#
}
\seealso{
  internal2tips wUniFrac
}
\keyword{internal}

