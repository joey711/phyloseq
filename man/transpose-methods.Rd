\docType{methods}
\name{t}
\alias{t}
\alias{t,otuTable-method}
\alias{t,phyloseq-method}
\title{Transpose \code{\link{otuTable-class}} or \code{\link{phyloseq-class}}}
\usage{
  t(x)
}
\arguments{
  \item{x}{An \code{otuTable} or
  \code{\link{phyloseq-class}}.}
}
\value{
  The class of the object returned by \code{t} matches the
  class of the argument, \code{x}. The \code{otuTable} is
  transposed, and \code{\link{speciesAreRows}} value is
  toggled.
}
\description{
  Extends the base transpose method, \code{\link[base]{t}}.
}
\examples{
data(GlobalPatterns)
otuTable(GlobalPatterns)
t( otuTable(GlobalPatterns) )
}

