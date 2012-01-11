\name{phyloseq-class}
\alias{phyloseq-class}
\title{The main experiment-level class for phyloseq data}
\description{
  Contains all component classes:
  \code{\link{otuTable-class}},
  \code{\link{sampleData-class}},
  \code{\link{taxonomyTable-class}} (\code{"taxTab"} slot),
  and \code{\link[ape]{phylo}}-class (\code{"tre"} slot).
  There are several advantages to storing your phylogenetic
  sequencing experiment as an instance of the phyloseq
  class, not the least of which is that it is easy to
  return to the data later and feel confident that the
  different data types ``belong'' to one another.
  Furthermore, the \code{\link{phyloseq}} constructor
  ensures that the different data components have
  compatible indices (e.g. species and samples), and
  performs the necessary trimming automatically when you
  create your ``experiment-level'' object. Downstream
  analyses are aware of which data classes they require --
  and where to find them -- often making your
  \code{phyloseq-class} object the only data argument to
  analysis and plotting functions (although there are many
  options and parameter arguments waiting for you).
}
\details{
  In the case of missing component data, the slots are set
  to \code{NULL}. As soon as a \code{phyloseq-class} object
  is to be updated with new component data (previously
  missing/\code{NULL} or not), the indices of all
  components are re-checked for compatibility and trimmed
  if necessary. This is to ensure by design that components
  describe the same taxa/samples, and also that these
  trimming/validity checks do not need to be repeated in
  downstream analyses.

  slots: \describe{ \item{otuTable}{a single object of
  class otuTable.} \item{samData}{ a single object of class
  sampleData.} \item{taxTab}{ a single object of class
  taxonomyTable.} \item{tre}{ a single object of class
  phylo, from the package ape} }
}
\seealso{
  The constructor, \code{\link{phyloseq}}, the merger
  \code{\link{merge_phyloseq}}, and also the component
  constructor/accessors \code{\link{otuTable}},
  \code{\link{sampleData}}, \code{\link{taxTab}}, and
  \code{\link{tre}}.
}

