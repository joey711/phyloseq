\docType{methods}
\name{taxa_names<-}
\alias{assign-taxa_names}
\alias{taxa_names<-}
\alias{taxa_names<-,ANY,ANY-method}
\alias{taxa_names<-,ANY,character-method}
\alias{taxa_names<-,otu_table,character-method}
\alias{taxa_names<-,phylo,character-method}
\alias{taxa_names<-,phyloseq,character-method}
\alias{taxa_names<-,taxonomyTable,character-method}
\alias{taxa_names<-,XStringSet,character-method}
\title{Replace OTU identifier names}
\usage{
  taxa_names(x) <- value
}
\arguments{
  \item{x}{(Required). An object defined by the
  \code{\link{phyloseq-package}} that describes OTUs in
  some way.}

  \item{value}{(Required). A character vector to replace
  the current \code{\link{taxa_names}}.}
}
\description{
  Replace OTU identifier names
}
\examples{
data("esophagus")
taxa_names(esophagus)
# plot_tree(esophagus, label.tips="taxa_names", ladderize="left")
taxa_names(esophagus) <- paste("OTU-", taxa_names(esophagus), sep="")
taxa_names(esophagus)
# plot_tree(esophagus, label.tips="taxa_names", ladderize="left")
## non-characters are first coerced to characters.
taxa_names(esophagus) <- 1:ntaxa(esophagus)
taxa_names(esophagus)
# plot_tree(esophagus, label.tips="taxa_names", ladderize="left")
## Cannot assign non-unique or differently-lengthed name vectors. Error.
# taxa_names(esophagus) <- sample(c(TRUE, FALSE), ntaxa(esophagus), TRUE)
# taxa_names(esophagus) <- sample(taxa_names(esophagus), ntaxa(esophagus)-5, FALSE)
}

