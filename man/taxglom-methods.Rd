\docType{methods}
\name{taxglom}
\alias{taxglom}
\alias{taxglom,otuTable,character-method}
\alias{taxglom,otuTable,taxonomyTable-method}
\alias{taxglom,phyloseq,ANY-method}
\title{Agglomerate taxa of the same type.}
\usage{
  taxglom(object, tax=NULL, taxlevel="Phylum")
}
\arguments{
  \item{object}{An \code{otuTable} or more complex object
  that contains an \code{otuTable}.}

  \item{tax}{Either a \code{taxonomyTable} object, or
  alternatively, a character vector specifying the desired
  taxonomic group of each taxa in \code{object}. If
  \code{tax} is a character vector, it must have length
  equal to the (original) number of taxa in \code{object},
  and each element must be named according to the taxa ID
  (that is, the result of \code{species.names(object)}). If
  \code{tax} is a character vector, than the
  \code{taxlevel} argument is ignored. If \code{object}
  already contains a \code{taxonomyTable} component in its
  \code{taxTab} slot, then the \code{tax} argument is
  ignored.}

  \item{taxlevel}{A single-element character specifying the
  taxonomic level (column name) in \code{tax}, the
  \code{taxonomyTable}, that you want to agglomerate over.
  The default value is \code{"Phylum"}. Note that this
  default may agglomerate too broadly for a given
  experiment, and the user is strongly encouraged to try
  different taxonomic levels.}
}
\value{
  A taxonomically-agglomerated object with class matching
  the class of \code{object}.
}
\description{
  This method merges species that are of the same
  taxaonomic level, using \code{mergespecies} similar in
  approach to \code{tipglom}, but using the categorical
  taxonomic levels rather than patristic distances on the
  phylogenetic tree.
}
\examples{
##
}

