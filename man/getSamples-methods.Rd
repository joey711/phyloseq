\docType{methods}
\name{getSamples}
\alias{getSamples}
\alias{getSamples,otuTable-method}
\alias{getSamples,phyloseq-method}
\title{Returns all abundance values for species \code{i}.}
\arguments{
  \item{physeq}{(Required). \code{\link{otuTable-class}},
  or \code{\link{phyloseq-class}}.}

  \item{i}{(Required). A single taxa/species/OTU ID for
  which you want to know the abundance in each sample.}
}
\value{
  An integer vector of the abundance values for each sample
  in \code{physeq} for species \code{i}
}
\description{
  This is a simple accessor function for investigating a
  single species-of-interest.
}
\examples{
data(esophagus)
species.names(esophagus)
getSamples(esophagus, "59_5_19")
}
\seealso{
  getSpecies species.names sample.names
}

