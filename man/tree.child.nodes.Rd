\name{tree.child.nodes}
\alias{tree.child.nodes}
\title{Return a list of a node's children.}
\usage{
  tree.child.nodes(phylo, node)
}
\arguments{
  \item{phylo,}{input phylo object}

  \item{node,}{integer index of the node to test}
}
\value{
  list, a list containing the integer indices of the nodes
  directly beneath the given node.
}
\description{
  Returns a list (not a vector!) of the node indices of the
  given node's direct children. Returns (-1, -1) if the
  given node is a leaf.
}
\author{
  Gregory Jordan \email{gjuggler@gmail.com}
}
\seealso{
  This code is borrowed directly, with permission, from the
  not-yet-officially-released package, \code{ggphylo},
  currently only available from GitHub at:
  \url{https://github.com/gjuggler/ggphylo}
}
\keyword{internal}

