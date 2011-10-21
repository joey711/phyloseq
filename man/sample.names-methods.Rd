\docType{methods}
\name{sample.names}
\alias{sample.names}
\alias{sample.names,otuTable-method}
\alias{sample.names,phyloseqFather-method}
\alias{sample.names,sampleMap-method}
\alias{sampleNames}
\title{Get the sample names of the samples described by an object.}
\usage{
  sample.names(x)
}
\arguments{
  \item{x}{(Required). An object among the set of classes
  defined by the phyloseq package that represent samples.}
}
\value{
  A character vector. The names of the samples in \code{x}.
}
\description{
  This method works on otuTable and sampleMap objects, as
  well as the more complex objects that represent samples
  in an experiment.
}
\examples{
#
# library("picante")
# data("phylocom")
# tree <- phylocom$phylo
# OTU1 <- otuTable(phylocom$sample, speciesAreRows=FALSE)
# sample.names(OTU1)
# ex1 <- phyloseq(OTU1, tree)
# sample.names(ex1)
}
\seealso{
  species.names nsamples
}

