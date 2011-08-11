\name{otuTree}
\alias{otutree}
\alias{otuTree}
\title{Define the otuTree class.}
\arguments{
  \item{object}{An object among the set of classes defined
  by the phyloseq package that contain both a phylogenetic
  tree and an otuTable.}
}
\value{
  A phyloseq-package object of class \code{otuTree}.
}
\description{
  Contains otuTable, phylo (tre) slots.

  Subset just the "otuTree" portion of a H.O. object.
}
\examples{
#
# library("picante")
# data("phylocom")
# tree <- phylocom$phylo
# OTU1 <- otuTable(phylocom$sample, speciesAreRows=FALSE)
# ex1 <- phyloseq(OTU1, tree)
# class(ex1)
# otuTree(ex1)
}
\seealso{
  phyloseq
}
\keyword{subset}

