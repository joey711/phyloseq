\name{plot_sample_network}
\alias{plot_sample_network}
\title{Plot sample-wise microbiome network (ggplot2)}
\usage{
  plot_sample_network(g, physeq=NULL, color=NULL,
  shape=NULL, point_size=4, alpha=1, label="value", hjust =
  1.35, line_weight=0.5, line_color=color, line_alpha=0.4,
  layout.method=layout.fruchterman.reingold)
}
\arguments{
  \item{g}{(Required). An
  \code{\link[igraph]{igraph}}-class object created either
  by the convenience wrapper
  \code{\link{make_sample_network}}, or directly by the
  tools in the igraph-package.}

  \item{physeq}{(Optional). Default \code{NULL}. A
  \code{\link{phyloseq-class}} object on which \code{g} is
  based.}

  \item{color}{(Optional). Default \code{NULL}. The name of
  the sample variable in \code{physeq} to use for color
  mapping of points (graph vertices).}

  \item{shape}{(Optional). Default \code{NULL}. The name of
  the sample variable in \code{physeq} to use for shape
  mapping. of points (graph vertices).}

  \item{point_size}{(Optional). Default \code{4}. The size
  of the vertex points.}

  \item{alpha}{(Optional). Default \code{1}. A value
  between 0 and 1 for the alpha transparency of the vertex
  points.}

  \item{label}{(Optional). Default \code{"value"}. The name
  of the sample variable in \code{physeq} to use for
  labelling the vertex points.}

  \item{hjust}{(Optional). Default \code{1.35}. The amount
  of horizontal justification to use for each label.}

  \item{line_weight}{(Optional). Default \code{0.3}. The
  line thickness to use to label graph edges.}

  \item{line_color}{(Optional). Default \code{color}. The
  name of the sample variable in \code{physeq} to use for
  color mapping of lines (graph edges).}

  \item{line_alpha}{(Optional). Default \code{0.4}. The
  transparency level for graph-edge lines.}

  \item{layout.method}{(Optional). Default
  \code{layout.fruchterman.reingold}. A function (closure)
  that determines the placement of the vertices for drawing
  a graph. Should be able to take an
  \code{\link{igraph}}-class as sole argument, and return a
  two-column coordinate matrix with \code{nrow} equal to
  the number of vertices. For possible options already
  included in \code{igraph}-package, see the others also
  described in the help file:

  \code{\link[igraph]{layout.fruchterman.reingold}}}
}
\value{
  A \code{\link{ggplot}}2 plot.
}
\description{
  A custom plotting function for displaying graph objects
  created by \code{\link[igraph]{igraph}} from a
  phylogenetic sequencing experiment
  (\code{\link{phyloseq-class}}), using advanced
  \code{\link[ggplot2]{ggplot}}2 formatting.
}
\examples{
data(enterotype)
ig <- make_sample_network(enterotype, FALSE, max.dist=0.3)
plot_sample_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.3, label=NULL)
# Change distance parameter
ig <- make_sample_network(enterotype, FALSE, max.dist=0.2)
plot_sample_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.3, label=NULL)
}
\references{
  Code modified from code now hosted on GitHub by Scott
  Chamberlain:
  \url{https://github.com/SChamberlain/gggraph}

  The code most directly used/modified was first posted
  here:
  \url{http://www.r-bloggers.com/basic-ggplot2-network-graphs/}
}
\seealso{
  \code{\link{make_sample_network}}
}

