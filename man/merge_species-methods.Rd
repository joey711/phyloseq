\docType{methods}
\name{merge_species}
\alias{merge_species}
\alias{merge_species,otuTable-method}
\alias{merge_species,phylo-method}
\alias{merge_species,phyloseq-method}
\alias{merge_species,sampleData-method}
\alias{merge_species,taxonomyTable-method}
\title{Merge a subset of the species in \code{x} into one species/taxa/OTU.}
\usage{
  merge_species(x, eqspecies, archetype=1)
}
\arguments{
  \item{x}{(Required). An object that describes species
  (taxa). This includes \code{\link{phyloseq-class}},
  \code{\link{otuTable-class}},
  \code{\link{taxonomyTable-class}},
  \code{\link[ape]{phylo}}.}

  \item{eqspecies}{(Required). The species names, or
  indices, that should be merged together. If
  \code{length(eqspecies) < 2}, then the object \code{x}
  will be returned safely unchanged.}

  \item{archetype}{The index of \code{eqspecies} indicating
  the species that should be kept (default is 1) to
  represent the summed/merged group of species/taxa/OTUs.
  If archetype is not an index or index-name in
  \code{eqspecies}, the first will be used, and the value
  in archetype will be used as the index-name for the new
  species.}
}
\value{
  The object, \code{x}, in its original class, but with the
  specified species merged into one entry in all relevant
  components.
}
\description{
  Takes as input an object that describes species/taxa
  (e.g. \code{\link{phyloseq-class}},
  \code{\link{otuTable-class}}, \code{\link{phylo-class}},
  \code{\link{taxonomyTable-class}}), as well as a vector
  of species that should be merged. It is intended to be
  able to operate at a low-level such that related methods,
  such as \code{\link{tipglom}} and \code{\link{taxglom}}
  can both reliably call \code{merge_species} for their
  respective purposes.
}
\examples{
#
# # data(phylocom)
# # tree <- phylocom$phylo
# # otu  <- otuTable(phylocom$sample, speciesAreRows=FALSE)
# # otutree0 <- phyloseq(otu, tree)
# # plot(otutree0)
# # otutree1 <- merge_species(otutree0, tree$tip.label[1:8], 2)
# # plot(otutree1)
}
\seealso{
  \code{\link{tipglom}}, \code{\link{taxglom}},
  \code{\link{merge_phyloseq}}, \code{\link{merge_samples}}
}

