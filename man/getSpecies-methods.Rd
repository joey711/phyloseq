\docType{methods}
\name{getSpecies}
\alias{getSpecies}
\alias{getSpecies,otuTable-method}
\alias{getSpecies,phyloseq-method}
\title{Returns all abundance values of sample \code{i}.}
\usage{
  getSpecies(physeq, i)
}
\arguments{
  \item{physeq}{(Required). \code{\link{otuTable-class}},
  or \code{\link{phyloseq-class}}.}

  \item{i}{(Required). A single sample for which you want
  to know the abundance of each species. Can be integer for
  index value, or sample name.}
}
\value{
  An integer vector of the abundance values for each
  species in \code{physeq} for sample \code{i}
}
\description{
  This is a simple accessor function for investigating a
  single sample-of-interest.
}
\examples{
data(esophagus)
sample.names(esophagus)
getSpecies(esophagus, "B")
}
\seealso{
  getSpecies species.names sample.names
}

