\docType{methods}
\name{t}
\alias{t}
\alias{t,otuTable-method}
\alias{t,phyloseqFather-method}
\title{Extension of the (base) transpose function for otuTable}
\usage{
  t(x)
}
\arguments{
  \item{x}{An \code{otuTable}, or a higher-order class that
  contains an \code{otuTable}.}
}
\value{
  The class of the object returned by \code{t} matches the
  class of the argument, \code{x}. The \code{otuTable} is
  now transposed, and \code{speciesAreRows} is toggled.
}
\description{
  Only the otuTable of a more complex object is transposed.
}
\examples{
#
}

