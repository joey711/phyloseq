\docType{methods}
\name{tipglom}
\alias{tipglom}
\alias{tipglom,phylo,ANY-method}
\alias{tipglom,phylo,otuTable-method}
\alias{tipglom,phyloseq,ANY-method}
\title{Agglomerate nearby tips in tree.}
\usage{
  tipglom(tree, OTU, speciationMinLength=0.02)
}
\arguments{
  \item{tree}{An object of class \code{phyloseq}, in which
  case the OTU argument can be ommitted. If, alternatively,
  \code{tree} is a \code{phylo} class, then \code{OTU} is
  required.}

  \item{OTU}{An otuTable object. Optional. Ignored if
  \code{tree} is a \code{phyloseq} object. If \code{tree}
  is a \code{phylo} object and \code{OTU} is provided, then
  return will be an \code{phyloseq} object.}

  \item{speciationMinLength}{The minimum branch length that
  separates taxa. All tips of the tree separated by a
  cophenetic distance smaller than
  \code{speciationMinLength} will be agglomerated. Default
  is 0.02}
}
\value{
  An object of class \code{phyloseq}. Output class matches
  the class of \code{tree}, unless it is a \code{phylo}
  object, in which case \code{tipglom} returns an
  \code{phyloseq} object.
}
\description{
  tipglom determines those groups of tips in the tree that
  are closer than the threshold provided as
  \code{speciationMinLength}. That is, all tips of the tree
  separated by a cophenetic distance smaller than
  \code{speciationMinLength} will be agglomerated using
  \code{mergespecies}. Can be used to create non-trivial
  OTU table, by agglomerating nearby tips.
}
\examples{
#
# # # data(phylocom)
# # # otu  <- otuTable(phylocom$sample, speciesAreRows=FALSE)
# # # x1   <- phyloseq(otu, phylocom$phylo)
# # # print(x1)
}

