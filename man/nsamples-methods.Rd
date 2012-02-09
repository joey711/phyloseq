\docType{methods}
\name{nsamples}
\alias{nsamples}
\alias{nsamples,ANY-method}
\alias{nsamples,otuTable-method}
\alias{nsamples,phyloseq-method}
\alias{nsamples,sampleData-method}
\title{Get the number of samples.}
\usage{
  nsamples(physeq)
}
\arguments{
  \item{physeq}{A \code{\link{phyloseq-class}},
  \code{\link{sampleData}}, or
  \code{\link{otuTable-class}}.}
}
\value{
  An integer indicating the total number of samples.
}
\description{
  Get the number of samples.
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
  \code{\link{species.names}}, \code{\link{sample.names}},
  \code{\link{nspecies}}
}

