\docType{methods}
\name{getSpecies}
\alias{getSpecies}
\alias{getSpecies,otuTable-method}
\alias{getSpecies,phyloseqFather-method}
\title{Returns the abundance values of sample \code{i} for
all species in \code{x}.}
\usage{
  getSpecies(x, i)
}
\arguments{
  \item{x}{otuTable, or H.O. object containing an
  otuTable.}

  \item{i}{a single sample for which you want to know the
  abundance of each species.}
}
\value{
  An integer vector of the abundance values for each
  species in \code{x} for sample \code{i}
}
\description{
  This is a simple accessor function for investigating a
  single sample-of-interest.
}
\examples{
#
}
\seealso{
  getSpecies species.names sample.names
}

