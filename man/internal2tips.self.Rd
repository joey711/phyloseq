\name{internal2tips.self}
\alias{internal2tips.self}
\title{Custom version of \code{\link[picante]{internal2tips}}}
\usage{
  internal2tips.self(phy, int.node, return.names = TRUE)
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
  Internal function for \code{\link{UniFrac}}.
}
\details{
  A modified version of the \code{\link{internal2tips}}
  function, such that when a tip is provided as
  \code{int.node}, that tip is returned. This is a more
  intuitive behavior than the original picante version,
  which returns NULL. This is currently used in
  \code{\link{UniFrac}}.
}
\seealso{
  \code{\link[picante]{internal2tips}}
  \code{\link{UniFrac}}
}
\keyword{internal}

