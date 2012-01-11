\docType{methods}
\name{makenetwork}
\alias{makenetwork}
\title{Create a taxa graph adjacency matrix from an \code{otuTable}.}
\usage{
  makenetwork(physeq, plotgraph=TRUE, community=TRUE,
    threshold=0, incommon=0.4, method="jaccard")
}
\arguments{
  \item{physeq}{(Required). An \code{\link{otuTable-class}}
  or \code{\link{phyloseq-class}} object.}

  \item{plotgraph}{A logical. Default TRUE. Indicates
  whether or not to plot the network}

  \item{community}{A logical. Default TRUE. If TRUE, RETURN
  will include the community groups.}

  \item{threshold}{A single-value positive integer. Default
  0. Indicates the number of individuals/observations
  required for presense to be counted.  For example,
  threshold > 1 means that species with 2 or more reads are
  considered present.}

  \item{incommon}{The fraction of co-occurring samples
  required for a pair of species to be given an edge.
  Default is 0.4. This parameter has a strong influence on
  the resulting network, and must be chosen carefully.}

  \item{method}{A character string indicating the default
  distance method parameter provided to
  \code{\link{vegdist}}. Default is ``jaccard''.}
}
\value{
  Creates a plot of the adjacency graph and optionally
  returns the community groups.
}
\description{
  By default, this uses a presence/absence criteria for
  taxa to determine if samples (vertices) are connected by
  edges, and then plot the result using the
  \code{\link[igraph]{igraph-package}}.
}
\examples{
## data(enterotype)
## makenetwork(enterotype)
}
\author{
  Susan Holmes
}

