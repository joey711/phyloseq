\docType{methods}
\name{tipglom}
\alias{tipglom}
\alias{tipglom,otuTree,ANY-method}
\alias{tipglom,phylo,ANY-method}
\alias{tipglom,phylo,otuTable-method}
\alias{tipglom,phylo4,ANY-method}
\title{Agglomerate nearby tips in tree.}
\arguments{
  \item{tree}{An object of class \code{otuTree} or its
  superclasses, in which case the OTU argument can be
  ommitted. If, alternatively, \code{tree} is a
  \code{phylo} or \code{phylo4} object, then \code{OTU} is
  required.}

  \item{OTU}{An otuTable object. Optional. Ignored if
  \code{tree} is an \code{otuTree} object or one of its
  superclasses. If \code{tree} is a \code{phylo} or
  \code{phylo4} object and \code{OTU} is provided, then
  return will be an \code{otuTree} object.}

  \item{speciationMinLength}{The minimum branch length that
  separates taxa. All tips of the tree separated by a
  cophenetic distance smaller than
  \code{speciationMinLength} will be agglomerated. Default
  is 0.02}
}
\value{
  An object of class \code{otuTree} or its superclasses.
  Output class matches the class of \code{tree}, unless it
  is a \code{phylo}/\code{phylo4} object, in which case
  \code{tipglom} returns an \code{otuTree} object.
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
# # # tree <- as(phylocom$phylo, "phylo4")
# # # x1   <- phyloseq(otu, tree)
# # # print(x1)
# # # library("phylobase")
# # # plot(tre(x1))
# # # x2 <- tipglom(x1, speciationMinLength=2.1)
# # # plot(tre(x2))
}

