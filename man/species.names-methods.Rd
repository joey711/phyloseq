\docType{methods}
\name{species.names}
\alias{species.names}
\alias{species.names,ANY-method}
\alias{species.names,otuTable-method}
\alias{species.names,phylo-method}
\alias{species.names,phyloseq-method}
\alias{species.names,sampleData-method}
\alias{species.names,taxonomyTable-method}
\title{Get species / taxa names.}
\usage{
  species.names(physeq)
}
\arguments{
  \item{physeq}{\code{\link{phyloseq-class}},
  \code{\link{otuTable-class}},
  \code{\link{taxonomyTable-class}}, or
  \code{\link[ape]{phylo}}}
}
\value{
  A character vector of the names of the species in
  \code{physeq}.
}
\description{
  Get species / taxa names.
}
\examples{
#
# # From "picante" package
# data("phylocom")
# tree <- phylocom$phylo
# OTU1 <- otuTable(phylocom$sample, speciesAreRows=FALSE)
# species.names(tree)
# species.names(OTU1)
# physeq1 <- phyloseq(OTU1, tree)
# species.names(physeq1)
}
\seealso{
  nspecies
}

