\docType{methods}
\name{phy_tree<-}
\alias{assign-phy_tree}
\alias{phy_tree<-}
\alias{phy_tree<-,phyloseq,phylo-method}
\alias{phy_tree<-,phyloseq,phyloseq-method}
\alias{tre<-}
\title{Assign a (new) phylogenetic tree to \code{x}}
\usage{
  phy_tree(x) <- value
}
\arguments{
  \item{x}{(Required). \code{\link{phyloseq-class}}}

  \item{value}{(Required). \code{\link{phylo-class}}, or
  \code{\link{phyloseq-class}}}
}
\description{
  Assign a (new) phylogenetic tree to \code{x}
}
\examples{
#
# data(GlobalPatterns)
# # An example of pruning to just the first 100 taxa in GlobalPatterns.
# ex2a <- prune_species(taxa_names(GlobalPatterns)[1:100], GlobalPatterns)
# # The following 3 lines produces an ex2b that is equal to ex2a
# ex2b <- GlobalPatterns
# tree <- prune_species(taxa_names(GlobalPatterns)[1:100], phy_tree(GlobalPatterns))
# phy_tree(ex2b) <- tree
# identical(ex2a, ex2b)
# print(ex2b)
# # Example adding a phylo tree from phyloseq class
# ex2c <- phyloseq(otu_table(ex2b), sample_data(ex2b), tax_table(ex2b))
# phy_tree(ex2c) <- ex2b
# identical(ex2b, ex2c)
}

