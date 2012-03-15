\docType{methods}
\name{tre}
\alias{tre}
\alias{tre,ANY-method}
\alias{tre,phylo-method}
\title{Get phylogenetic tree from object.}
\usage{
  tre(physeq, errorIfNULL=TRUE)
}
\arguments{
  \item{physeq}{(Required). An instance of phyloseq-class
  that contains a phylogenetic tree. If physeq is a
  phylogenetic tree (a component data class), then it is
  returned as-is.}

  \item{errorIfNULL}{(Optional). Logical. Should the
  accessor stop with an error if the slot is empty
  (\code{NULL})? Default \code{TRUE}.}
}
\value{
  The \code{\link[ape]{phylo}}-class object contained
  within \code{physeq}; or NULL if \code{physeq} does not
  have a tree. This method stops with an error in the
  latter NULL case be default, which can be over-ridden by
  changing the value of \code{errorIfNULL} to \code{FALSE}.
}
\description{
  This is the main method suggested for accessing the
  phylogenetic tree, (\code{\link[ape]{phylo}}-class) from
  a \code{\link{phyloseq-class}}. Like other accessors (see
  See Also, below), the default behavior of this method is
  to stop with an error if \code{physeq} is a
  \code{phyloseq-class} but does not contain a phylogenetic
  tree.
}
\details{
  Note that the tip labels should be named to match the
  \code{species.names} of the other objects to which it is
  going to be paired. The \code{\link{phyloseq}}
  constructor automatically checks for exact agreement in
  the set of species described by the phlyogenetic tree and
  the other components (taxonomyTable, otuTable), and trims
  as-needed. Thus, the tip.labels in a phylo object must be
  named to match the results of \code{\link{species.names}}
  of the other objects to which it will ultimately be
  paired.
}
\examples{
# data(GlobalPatterns)
# tre(GlobalPatterns)
}
\seealso{
  \code{\link{otuTable}}, \code{\link{sampleData}},
  \code{\link{taxTab}} \code{\link{phyloseq}},
  \code{\link{merge_phyloseq}}
}

