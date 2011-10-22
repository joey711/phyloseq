\docType{methods}
\name{wUniFrac}
\alias{wUniFrac}
\alias{wUniFrac,otuTable,phylo-method}
\alias{wUniFrac,otuTree,ANY-method}
\title{Calculate weighted UniFrac for all samples in an OTU table.}
\usage{
  wUniFrac(OTU, tree, normalized=TRUE, num_cores=1)
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

  \item{num_cores}{(Optional). Integer. The number of CPU
  cores to use in calculation. This is not checked for
  validity. User must specify a valid number. If value is
  \code{<=1}, then a serial implementation is used. Default
  is 1.}
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

