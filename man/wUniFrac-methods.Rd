\docType{methods}
\name{wUniFrac}
\alias{wUniFrac}
\alias{wUniFrac,otuTable,phylo-method}
\alias{wUniFrac,otuTree,ANY-method}
\title{Calculate weighted UniFrac for all samples in an OTU table.}
\arguments{
  \item{OTU}{otuTable, or a more-complex object that
  contains an otuTable.}

  \item{tree}{object of class \code{phylo}}

  \item{normalized}{Logical. Should the output be
  normalized such that values range from 0 to 1 independent
  of branch length values? Default is \code{TRUE}.}
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

  Calculate weighted UniFrac from an abundance matrix and a
  tree.

  Calculate weighted UniFrac from an abundance matrix and a
  tree.
}
\examples{
#
}
\seealso{
  UniFrac unifrac vegdist
}

