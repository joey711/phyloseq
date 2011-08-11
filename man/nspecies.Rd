\name{nspecies}
\alias{nspecies}
\title{Get the number of taxa/species in an object.}
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
# library("picante")
# data("phylocom")
# tree <- phylocom$phylo
# nspecies(tree)
}
\seealso{
  species.names
}
\keyword{OTU}
\keyword{richness}
\keyword{species}
\keyword{taxa}

