\docType{methods}
\name{tip_glom}
\alias{tip_glom}
\alias{tip_glom,phylo,ANY-method}
\alias{tip_glom,phylo,otu_table-method}
\alias{tip_glom,phyloseq,ANY-method}
\alias{tipglom}
\title{Agglomerate closely-related taxa using single-linkage clustering.}
\usage{
  tip_glom(tree, OTU, speciationMinLength=0.02)

  tipglom(tree, OTU, speciationMinLength = 0.02)
}
\arguments{
  \item{tree}{\code{\link{phyloseq-class}}, containing an
  OTU Table and phylogenetic tree. If, alternatively,
  \code{tree} is a \code{\link{phylo-class}}, then
  \code{OTU} is required.}

  \item{OTU}{An otu_table object. Optional. Ignored if
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
  object, in which case \code{tip_glom} returns an
  \code{phyloseq} object.
}
\description{
  All tips of the tree separated by a cophenetic distance
  smaller than \code{speciationMinLength} will be
  agglomerated into one taxa using \code{merge_taxa}.
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
data(esophagus); ntaxa(esophagus); plot(phy_tree(esophagus))
phy_tree(esophagus)$edge.length
x3 <- tip_glom(esophagus, speciationMinLength = 0.20)
ntaxa(x3); plot(phy_tree(x3))
}

