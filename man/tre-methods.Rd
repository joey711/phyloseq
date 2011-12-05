\docType{methods}
\name{tre}
\alias{tre}
\alias{tre,ANY-method}
\alias{tre,phylo-method}
\title{Access tre slot, or check/coerce to phylo4 class.}
\usage{
  tre(object, errorIfNULL=TRUE)
}
\arguments{
  \item{object}{(Required). An instance of phyloseq-class
  that contains a phylogenetic tree. If object is a
  phylogenetic tree (a component data class), then it is
  returned as-is.}

  \item{errorIfNULL}{(Optional). Logical. Should the
  accessor stop with an error if the slot is empty
  (\code{NULL})? Default \code{TRUE}.}
}
\value{
  A phylo4 object. It is grabbed from the tre-slot. If
  object is not a phylogenetic tree and does not contain a
  tre-slot, then NULL is returned. This is a convenience
  wrapper of the \code{\link{access}} function.
}
\description{
  \code{tre()} is an accessor OR coercion method. This is
  the main method suggested for accessing/subsetting the
  phylogenetic tree (phylo4 class) from a more complex
  object. \code{tre} is the slot name for trees. that holds
  the phylo4-class object in a multi-component
  phyloseq-package object.
}
\details{
  Note that the tip labels should be named to match the
  \code{species.names} of the other objects to which it is
  going to be paired. The initialization methods for more
  complex objects automatically check for exact agreement
  in the set of species described by the phlyogenetic tree
  and the other components (taxonomyTable, otuTable). They
  also trim accordingly. Thus, the tip.labels in a phylo
  object must be named to match the \code{species.names} of
  the other objects to which it will ultimately be paired.
}
\seealso{
  otuTable sampleMap taxTab phyloseq
}

