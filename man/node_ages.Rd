\name{node_ages}
\alias{node_ages}
\title{Return the ages (length to root node) of each node in a phylogenetic tree}
\usage{
  node_ages(phy)
}
\arguments{
  \item{phy}{(Required). A phylogenetic tree in
  \code{\link[ape]{phylo}} format.}
}
\value{
  A numeric vector
}
\description{
  This is a reimplementation of the function
  \code{node.age} from
  \href{http://cran.at.r-project.org/web/packages/picante/index.html}{the
  picante package}. In testing, the more-vectorized version
  of the computation included here is roughly 3-fold faster
  than \code{node.age} and with identical results. Unlike
  the picante-package implementation, this function does
  not return a modified version of the phylogenetic tree
  input, rather it returns only the newly-computed vector
  of ages for each node.
}
\keyword{internal}

