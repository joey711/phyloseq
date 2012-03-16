\docType{methods}
\name{tipglom}
\alias{tipglom}
\alias{tipglom,phylo,ANY-method}
\alias{tipglom,phylo,otuTable-method}
\alias{tipglom,phyloseq,ANY-method}
\title{Agglomerate closely-related taxa using single-linkage clustering.}
\usage{
  tipglom(tree, OTU, speciationMinLength=0.02)
}
\arguments{
  \item{tree}{\code{\link{phyloseq-class}}, containing an
  OTU Table and phylogenetic tree. If, alternatively,
  \code{tree} is a \code{\link{phylo-class}}, then
  \code{OTU} is required.}

  \item{OTU}{An otuTable object. Optional. Ignored if
  \code{tree} is a \code{\link{phyloseq-class}} object. If
  \code{tree} is a \code{phylo} object and \code{OTU} is
  provided, then return will be an \code{phyloseq} object.}

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
  All tips of the tree separated by a cophenetic distance
  smaller than \code{speciationMinLength} will be
  agglomerated into one taxa using \code{merge_species}.
}
\details{
  Can be used to create a non-trivial OTU Table, if a
  phylogenetic tree is available.

  For now, a simple, ``greedy'', single-linkage clustering
  is used. In future releases it should be possible to
  specify different clustering approaches available in
  \code{R}, in particular, complete-linkage clustering
  appears to be used more commonly for OTU clustering
  applications.
}
\examples{
#
# # # data(phylocom)
# # # otu  <- otuTable(phylocom$sample, speciesAreRows=FALSE)
# # # x1   <- phyloseq(otu, phylocom$phylo)
# # # print(x1); par(mfrow=c(2, 1)); plot(tre(x1))
# # # x2 <- tipglom(x1, speciationMinLength = 2.5)
# # # plot(tre(x2))
# # # ## Try on example datset 1
# # # data(GlobalPatterns); nspecies(GlobalPatterns)
# # # ex7 <- tipglom(GlobalPatterns, speciationMinLength = 0.05)
# # # nspecies(ex7)
# data(esophagus); nspecies(esophagus); par(mfrow=c(2, 1)); plot(tre(esophagus))
# tre(esophagus)$edge.length
# x3 <- tipglom(esophagus, speciationMinLength = 0.20)
# nspecies(x3); plot(tre(x3))
}

