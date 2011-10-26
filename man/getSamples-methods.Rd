\docType{methods}
\name{getSamples}
\alias{getSamples}
\alias{getSamples,otuTable-method}
\alias{getSamples,phyloseqFather-method}
\title{Returns the abundance values of species \code{i} for
all samples in \code{x}.}
\arguments{
  \item{x}{otuTable, or H.O. object containing an
  otuTable.}

  \item{i}{a single taxa/species/OTU ID for which you want
  to know the abundance in each sample.}
}
\value{
  An integer vector of the abundance values for each sample
  in \code{x} for species \code{i}
}
\description{
  This is a simple accessor function for investigating a
  single species-of-interest.
}
\examples{
#
}
\seealso{
  getSpecies species.names sample.names
}

