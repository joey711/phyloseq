\name{internal2tips_self}
\alias{internal2tips_self}
\title{Custom version of picante's \code{internal2tips}}
\usage{
  internal2tips_self(phy, int.node, return.names = TRUE)
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
  A modified version of the \code{internal2tips} function
  from
  \href{http://cran.at.r-project.org/web/packages/picante/index.html}{the
  picante package}, such that when a tip is provided as
  \code{int.node}, that tip is returned. This is a more
  intuitive behavior than the original picante version,
  which returns \code{NULL}. This is currently used in
  \code{\link{UniFrac}}.
}
\seealso{
  \code{internal2tips}

  \code{\link{UniFrac}}
}
\keyword{internal}

