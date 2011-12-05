\docType{methods}
\name{nsamples}
\alias{nsamples}
\alias{nsamples,ANY-method}
\alias{nsamples,phyloseq-method}
\title{Get the number of samples described by an object.}
\usage{
  nsamples(object)
}
\arguments{
  \item{object}{An object among the set of classes defined
  by the phyloseq package that represent samples.}
}
\value{
  An integer indicating the total number of samples.
}
\description{
  This method works on otuTable and sampleMap objects, as
  well as the more complex objects that represent samples
  in an experiment.
}
\examples{
#
# # From "picante" package
# data("phylocom")
# tree <- phylocom$phylo
# OTU1 <- otuTable(phylocom$sample, speciesAreRows=FALSE)
# nsamples(OTU1)
# physeq1 <- phyloseq(OTU1, tree)
# nsamples(physeq1)
}
\seealso{
  nspecies sample.names
}

