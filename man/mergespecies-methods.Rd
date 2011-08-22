\docType{methods}
\name{mergespecies}
\alias{"mergespecies,otuTable-method"}
\alias{"mergespecies,otuTree-method"}
\alias{"mergespecies,otuTree4-method"}
\alias{"mergespecies,phylo-method"}
\alias{"mergespecies,phylo4-method"}
\alias{"mergespecies,phyloseq-method"}
\alias{"mergespecies,phyloseqTaxTree-method"}
\alias{"mergespecies,phyloseqTree-method"}
\alias{mergespecies}
\title{Merge a subset of the species in \code{x} into one species/taxa/OTU.}
\arguments{
  \item{x}{The otuTable or phylo object, or a higher-order
  object that contains an otuTable or tree.}

  \item{eqspecies}{The species names, or indices, that
  should be merged together. If \code{length(eqspecies) <
  2}, then the object \code{x} will be returned safely
  unchanged.}

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
  \code{mergespecies} is a method that takes as input an
  otuTable (or higher object) and a vector of species that
  should be merged. It is intended to be able to operate at
  a low-level such that related methods, such as tipglom
  and taxglom can both reliably call mergespecies for their
  respective purposes.

  Merge a subset of the taxa in an otuTable.

  Merge a subset of the taxa in a phylo-class tree object.

  Merge a subset of the taxa in a phylo4-class tree object.

  Note to developers: The class name otuTree4 is a
  placeholder. plan would be to migrate all ape-derived
  "phylo" classes to \code{phylo4} of phylobase package
  (which has a namespace). Lots more potential with phylo4,
  phylo4d, etc. E.g. merge could store the "lost" tips as
  data (phylo4d) associated with each tip. More effective
  than making separate list object to keep track.

  Merge a subset of the taxa in an otuTree object.

  Merge a subset of the taxa in an phyloseq object.

  Merge a subset of the taxa in an phyloseqTree object.

  Merge a subset of the taxa in an phyloseqTaxTree object.
}
\examples{
#
## data(ex1)
## makenetwork(otuTable(ex1), TRUE)
}
\seealso{
  tipglom taxglom merge_phyloseq
}

