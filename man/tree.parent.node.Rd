\name{tree.parent.node}
\alias{tree.parent.node}
\title{Extracts the parent node index for the given node.}
\usage{
  tree.parent.node(phylo, node)
}
\arguments{
  \item{phylo,}{input phylo object}

  \item{node,}{integer index of the node whose parent is
  desired}
}
\value{
  integer, the index of the parent node or -1 if the given
  node is root.
}
\description{
  Returns -1 if the node is root. Return the index of the
  node directly above the given node. Returns -1 if the
  given node is root.
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

