\name{tree.branch.length}
\alias{tree.branch.length}
\title{Extracts the length of the branch above the given node.}
\usage{
  tree.branch.length(phylo, node)
}
\arguments{
  \item{phylo}{input phylo object}

  \item{node}{integer, the node's index}
}
\value{
  numeric, the branch length of the edge leading to the
  given node. May be NA.
}
\description{
  Returns 0 if the node is root.
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

