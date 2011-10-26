\docType{methods}
\name{makenetwork}
\alias{makenetwork}
\alias{makenetwork,otuTable-method}
\alias{makenetwork,phyloseqFather-method}
\title{Create a taxa graph adjacency matrix from an \code{otuTable}.}
\arguments{
  \item{abund}{The species abundance table (otuTable) or a
  more complex phyloseq-package object that contains an
  otuTable.}

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
  Function that inputs an abundance table as a matrix and
  makes a graph adjacency matrix ready to plot with the
  igraph package.
}
\examples{
#
## data(ex1)
## makenetwork(otuTable(ex1), TRUE)
}
\author{
  Susan Holmes
}

