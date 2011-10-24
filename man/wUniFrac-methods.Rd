\docType{methods}
\name{wUniFrac}
\alias{wUniFrac}
\alias{wUniFrac,otuTable,phylo-method}
\alias{wUniFrac,otuTree,ANY-method}
\title{Calculate weighted UniFrac for all samples in an OTU table.}
\usage{
  wUniFrac(OTU, tree, normalized=TRUE, parallel=FALSE)
}
\arguments{
  \item{OTU}{(Required). \code{otuTable}, or an
  \code{otuTree}. If you have instead a simple matrix of
  abundances, see \code{\link{otuTable}} for coercing it to
  the \code{otuTable} class.}

  \item{tree}{(Optional). Object of class \code{phylo}. Not
  necessary (and ignored) if \code{OTU} is an object that
  also contains a tree (\code{otuTree} class, or its
  children).}

  \item{normalized}{(Optional). Logical. Should the output
  be normalized such that values range from 0 to 1
  independent of branch length values? Default is
  \code{TRUE}.}

  \item{parallel}{(Optional). Logical. Should execute
  calculation in parallel, using multiple CPU cores
  simultaneously? This can dramatically hasten the
  computation time for this function. However, it also
  requires that the user has registered a parallel
  ``backend'' prior to calling this function. Default is
  \code{FALSE}. If FALSE, wUniFrac will register a serial
  backend so that \code{foreach::\%dopar\%} does not throw
  a warning.}
}
\value{
  a sample-by-sample distance matrix, suitable for NMDS,
  etc.
}
\description{
  This function calculates the weighted-UniFrac distance
  for all sample-pairs in a species-abundance table using
  the abundances and a phylogenetic tree. If \code{OTU} is
  a more complex object that already contains a
  phylogenetic tree and abundance table, then the argument
  \code{tree} is not necessary and will be ignored.
}
\examples{
#
}
\seealso{
  UniFrac unifrac vegdist
}

