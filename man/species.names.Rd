\name{species.names}
\alias{species.names}
\title{Get the species / taxa names from an object.}
\arguments{
  \item{object}{An object among the set of classes defined
  by the phyloseq package that represent species.}
}
\value{
  A character vector of the names of the species in
  \code{object}.
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
# OTU1 <- otuTable(phylocom$sample, speciesAreRows=FALSE)
# species.names(tree)
# species.names(OTU1)
# physeq1 <- phyloseq(OTU1, tree)
# species.names(physeq1)
}
\seealso{
  nspecies
}
\keyword{OTU}
\keyword{richness}
\keyword{species}
\keyword{taxa}

