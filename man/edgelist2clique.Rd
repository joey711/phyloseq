\name{edgelist2clique}
\alias{edgelist2clique}
\title{Convert edgelist hash-table to clique list}
\usage{
  edgelist2clique(EdgeList)
}
\arguments{
  \item{EdgeList}{a 2-column table of vertices where each
  row represents an edge.}
}
\value{
  A list, where each element is a character vector of tips
  that should are in the same clique.
}
\description{
  Agglomerating function to convert an edgelist -- which is
  a 2-column table of vertices where each row represents an
  edge -- to a list of cliques, in which each clique is
  represented by a character vector of the vertex labels of
  the vertices that are members of the clique. This
  algorithm is perfectly greedy, such that the only
  requirement for inclusion in a clique is an edge to any
  of the other members of that clique.
}
\examples{
#
# edgelist2clique(get.edgelist(ig))
}
\keyword{internal}

