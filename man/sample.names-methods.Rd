\docType{methods}
\name{sample.names}
\alias{sample.names}
\alias{sample.names,ANY-method}
\alias{sample.names,otuTable-method}
\alias{sample.names,phyloseq-method}
\alias{sample.names,sampleData-method}
\alias{sampleNames}
\title{Get sample names.}
\arguments{
  \item{physeq}{(Required). A \code{\link{phyloseq-class}},
  \code{\link{sampleData}}, or
  \code{\link{otuTable-class}}.}
}
\value{
  A character vector. The names of the samples in
  \code{physeq}.
}
\description{
  Get sample names.
}
\examples{
#
# # From "picante" package
# data("phylocom")
# tree <- phylocom$phylo
# OTU1 <- otuTable(phylocom$sample, speciesAreRows=FALSE)
# sample.names(OTU1)
# ex1 <- phyloseq(OTU1, tree)
# sample.names(ex1)
}
\seealso{
  \code{\link{species.names}}, \code{\link{nsamples}}
}

