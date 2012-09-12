\docType{methods}
\name{ntaxa}
\alias{nspecies}
\alias{ntaxa}
\alias{ntaxa,ANY-method}
\alias{ntaxa,otu_table-method}
\alias{ntaxa,phylo-method}
\alias{ntaxa,phyloseq-method}
\alias{ntaxa,taxonomyTable-method}
\title{Get the number of taxa/species.}
\usage{
  ntaxa(physeq)

  nspecies(physeq)
}
\arguments{
  \item{physeq}{\code{\link{phyloseq-class}},
  \code{\link{otu_table-class}},
  \code{\link{taxonomyTable-class}}, or
  \code{\link[ape]{phylo}}}
}
\value{
  An integer indicating the number of taxa / species.
}
\description{
  Get the number of taxa/species.
}
\examples{
#
# # From "picante" package
# data("phylocom")
# tree <- phylocom$phylo
# ntaxa(tree)
}
\seealso{
  taxa_names
}

