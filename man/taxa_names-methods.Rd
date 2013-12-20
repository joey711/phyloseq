\docType{methods}
\name{taxa_names}
\alias{taxa_names}
\alias{taxa_names,ANY-method}
\alias{taxa_names,otu_table-method}
\alias{taxa_names,phylo-method}
\alias{taxa_names,phyloseq-method}
\alias{taxa_names,sample_data-method}
\alias{taxa_names,taxonomyTable-method}
\alias{taxa_names,XStringSet-method}
\title{Get species / taxa names.}
\usage{
  taxa_names(physeq)
}
\arguments{
  \item{physeq}{\code{\link{phyloseq-class}},
  \code{\link{otu_table-class}},
  \code{\link{taxonomyTable-class}}, or
  \code{\link[ape]{phylo}}}
}
\value{
  A character vector of the names of the species in
  \code{physeq}.
}
\description{
  Get species / taxa names.
}
\examples{
#
data("esophagus")
tree <- phy_tree(esophagus)
OTU1 <- otu_table(esophagus)
taxa_names(tree)
taxa_names(OTU1)
physeq1 <- phyloseq(OTU1, tree)
taxa_names(physeq1)
}
\seealso{
  ntaxa
}

