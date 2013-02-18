\docType{methods}
\name{merge_taxa}
\alias{merge_species}
\alias{merge_taxa}
\alias{merge_taxa,otu_table-method}
\alias{merge_taxa,phylo-method}
\alias{merge_taxa,phyloseq-method}
\alias{merge_taxa,sample_data-method}
\alias{merge_taxa,taxonomyTable-method}
\alias{merge_taxa,XStringSet-method}
\title{Merge a subset of the species in \code{x} into one species/taxa/OTU.}
\usage{
  merge_taxa(x, eqtaxa, archetype=1)

  merge_species(x, eqtaxa, archetype = 1L)
}
\arguments{
  \item{x}{(Required). An object that describes species
  (taxa). This includes \code{\link{phyloseq-class}},
  \code{\link{otu_table-class}},
  \code{\link{taxonomyTable-class}},
  \code{\link[ape]{phylo}}.}

  \item{eqtaxa}{(Required). The species names, or indices,
  that should be merged together. If \code{length(eqtaxa) <
  2}, then the object \code{x} will be returned safely
  unchanged.}

  \item{archetype}{The index of \code{eqtaxa} indicating
  the species that should be kept (default is 1) to
  represent the summed/merged group of species/taxa/OTUs.
  If archetype is not an index or index-name in
  \code{eqtaxa}, the first will be used, and the value in
  archetype will be used as the index-name for the new
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
  \code{\link{otu_table-class}}, \code{\link{phylo-class}},
  \code{\link{taxonomyTable-class}}), as well as a vector
  of species that should be merged. It is intended to be
  able to operate at a low-level such that related methods,
  such as \code{\link{tip_glom}} and \code{\link{tax_glom}}
  can both reliably call \code{merge_taxa} for their
  respective purposes.
}
\examples{
#
data(phylocom)
tree <- phylocom$phylo
otu  <- otu_table(phylocom$sample, taxa_are_rows=FALSE)
otutree0 <- phyloseq(otu, tree)
plot_tree(otutree0)
otutree1 <- merge_taxa(otutree0, 1:8, 2)
plot_tree(otutree1)
}
\seealso{
  \code{\link{tip_glom}}, \code{\link{tax_glom}},
  \code{\link{merge_phyloseq}}, \code{\link{merge_samples}}
}

