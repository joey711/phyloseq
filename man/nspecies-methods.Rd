\docType{methods}
\name{nspecies}
\alias{nspecies}
\alias{nspecies,otuTable-method}
\alias{nspecies,phylo-method}
\alias{nspecies,phylo4-method}
\alias{nspecies,phyloseq-method}
\alias{nspecies,taxonomyTable-method}
\title{Get the number of taxa/species in an object.}
\usage{
  nspecies(object)
}
\arguments{
  \item{object}{An object among the set of classes defined
  by the phyloseq package that represent species.}
}
\value{
  An integer indicating the number of taxa / species.
}
\description{
  This method works on otuTable, taxonomyTable objects, or
  the more complex objects that represent species, as well
  as phylogenetic trees ``phylo''.
}
\examples{
#
# # From "picante" package
# data("phylocom")
# tree <- phylocom$phylo
# nspecies(tree)
}
\seealso{
  species.names
}

