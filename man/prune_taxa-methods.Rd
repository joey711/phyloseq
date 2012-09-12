\name{prune_taxa}
\alias{prune_species}
\alias{prune_taxa}
\alias{prune_taxa,NULL,ANY-method}
\alias{prune_taxa,character,otu_table-method}
\alias{prune_taxa,character,phylo-method}
\alias{prune_taxa,character,phyloseq-method}
\alias{prune_taxa,character,sample_data-method}
\alias{prune_taxa,character,taxonomyTable-method}
\alias{prune_taxa,logical,ANY-method}
\title{Prune unwanted OTUs / taxa from a phylogenetic object.}
\usage{
  prune_species(taxa, x)
}
\arguments{
  \item{taxa}{(Required). A character vector of the taxa in
  object x that you want to keep -- OR alternatively -- a
  logical vector where the kept taxa are TRUE, and length
  is equal to the number of taxa in object x. If
  \code{taxa} is a named logical, the taxa retained are
  based on those names. Make sure they are compatible with
  the \code{taxa_names} of the object you are modifying
  (\code{x}).}

  \item{x}{(Required). A phylogenetic object, including
  \code{phylo} trees, as well as all phyloseq classes that
  represent taxa. If the function \code{\link{taxa_names}}
  returns a non-\code{NULL} value, then your object can be
  pruned by this function.}
}
\value{
  The class of the object returned by \code{prune_taxa}
  matches the class of the argument, \code{x}.
}
\description{
  An S4 Generic method for removing (pruning) unwanted
  OTUs/taxa from phylogenetic objects, including
  phylo-class trees, as well as native phyloseq package
  objects. This is particularly useful for pruning a
  phyloseq object that has more than one component that
  describes OTUs. Credit: the \code{phylo}-class version is
  adapted from \code{\link[picante]{prune.sample}}.
}
\examples{
#
## testOTU <- otu_table(matrix(sample(1:50, 25, replace=TRUE), 5, 5), taxa_are_rows=FALSE)
## f1  <- filterfun_sample(topk(2))
## wh1 <- genefilter_sample(testOTU, f1, A=2)
## wh2 <- c(T, T, T, F, F)
## prune_taxa(wh1, testOTU)
## prune_taxa(wh2, testOTU)
##
## tax_table1 <- tax_table(matrix("abc", 5, 5))
## prune_taxa(wh1, tax_table1)
## prune_taxa(wh2, tax_table1)
}
\seealso{
  \code{\link{prune_taxa}}

  \code{\link[picante]{prune.sample}}
}

