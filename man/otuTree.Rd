\name{otuTree}
\alias{otutree}
\alias{otuTree}
\title{Create a subset object with just the \code{otuTable} and \code{tre} slots.}
\usage{
  otuTree(object)
}
\arguments{
  \item{object}{An object among the set of classes defined
  by the phyloseq package that contain both a phylogenetic
  tree and an otuTable.}
}
\value{
  A phyloseq-package object of class \code{otuTree}.
}
\description{
  Create a subset object with just the \code{otuTable} and
  \code{tre} slots.
}
\examples{
#
# # From "picante" package
# data("phylocom")
# tree <- phylocom$phylo
# OTU1 <- otuTable(phylocom$sample, speciesAreRows=FALSE)
# ex1 <- phyloseq(OTU1, tree)
# class(ex1)
# otuTree(ex1)
}

