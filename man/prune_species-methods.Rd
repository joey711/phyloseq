\name{prune_species}
\alias{prune_species}
\alias{prune_species,character,otuTable-method}
\alias{prune_species,character,phylo-method}
\alias{prune_species,character,phyloseq-method}
\alias{prune_species,character,sampleData-method}
\alias{prune_species,character,taxonomyTable-method}
\alias{prune_species,logical,ANY-method}
\alias{prune_species,NULL,ANY-method}
\title{Prune unwanted species / taxa from a phylogenetic object.}
\arguments{
  \item{species}{(Required). A character vector of the
  species in object x that you want to keep -- OR
  alternatively -- a logical vector where the kept species
  are TRUE, and length is equal to the number of species in
  object x. If \code{species} is a named logical, the
  species retained is based on those names. Make sure they
  are compatible with the \code{species.names} of the
  object you are modifying (\code{x}).}

  \item{x}{(Required). A phylogenetic object, including
  \code{phylo} trees, as well as all phyloseq classes that
  represent taxa / species. If the function
  \code{\link{species.names}} returns a non-\code{NULL}
  value, then your object can be pruned by this function.}
}
\value{
  The class of the object returned by \code{prune_species}
  matches the class of the argument, \code{x}.
}
\description{
  An S4 Generic method for removing (pruning) unwanted taxa
  from phylogenetic objects, including phylo-class trees,
  as well as native phyloseq package objects. This is
  particularly useful for pruning a phyloseq object that
  has more than one component that describes species. The
  \code{phylo}-class version is adapted from
  \code{picante::prune.samples}.
}
\examples{
#
## testOTU <- otuTable(matrix(sample(1:50, 25, replace=TRUE), 5, 5), speciesAreRows=FALSE)
## f1  <- filterfunSample(topk(2))
## wh1 <- genefilterSample(testOTU, f1, A=2)
## wh2 <- c(T, T, T, F, F)
## prune_species(wh1, testOTU)
## prune_species(wh2, testOTU)
##
## taxtab1 <- taxTab(matrix("abc", 5, 5))
## prune_species(wh1, taxtab1)
## prune_species(wh2, taxtab1)
}

