\name{tree.get.tags}
\alias{tree.get.tags}
\title{Retrieves a list of all tags for the given node.}
\usage{
  tree.get.tags(phylo, node)
}
\arguments{
  \item{phylo}{input phylo object}

  \item{node}{the node index for the desired tags}
}
\value{
  list containing all tags associated with this node, if
  tags exist; empty list otherwise.
}
\description{
  Retrieves a list of all tags for the given node.
}
\examples{
# tree <- tree.read('((a,b[&&NHX:foo=bar]),c);')
# tree.get.tags(tree, tree.node.with.label(tree, 'b')) # foo => bar
#
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

