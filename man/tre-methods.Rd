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
  The \code{\link{phylo4-class}} object contained within
  \code{object}; or NULL if \code{object} does not have a
  tree. This method stops with an error in the latter NULL
  case be default, which can be over-ridden by changing the
  value of \code{errorIfNULL} to \code{FALSE}.
}
\description{
  \code{tre()} is an accessor OR coercion method. This is
  the main method suggested for accessing the phylogenetic
  tree (\code{\link{phylo4-class}}) from a
  \code{\link{phyloseq-class}}. Like other accessors (see
  See Also, below), the default behavior of this method is
  to stop with an error if \code{object} is a
  \code{phyloseq-class} but does not contain a phylogenetic
  tree.
}
\details{
  Note that the tip labels should be named to match the
  \code{species.names} of the other objects to which it is
  going to be paired. The initialization methods for more
  complex objects automatically check for exact agreement
  in the set of species described by the phlyogenetic tree
  and the other components (taxonomyTable, otuTable). They
  also trim accordingly. Thus, the tip.labels in a phylo
  object must be named to match the results of
  \code{\link{species.names}} of the other objects to which
  it will ultimately be paired.
}
\examples{
# data(ex1)
# tre(ex1)
}
\seealso{
  \code{\link{otuTable}}, \code{\link{sampleMap}},
  \code{\link{taxTab}} \code{\link{phyloseq}},
  \code{\link{merge_phyloseq}}
}

