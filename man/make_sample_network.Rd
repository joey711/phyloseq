\name{make_sample_network}
\alias{make_sample_network}
\title{Make sample-wise microbiome network (igraph)}
\usage{
  make_sample_network(physeq, dist.fun="jaccard", max.dist
  = 0.4, keep.isolates=FALSE, ...)
}
\arguments{
  \item{physeq}{(Required). Default \code{NULL}. A
  \code{\link{phyloseq-class}} object, or
  \code{\link{otuTable-class}} object, on which \code{g} is
  based. \code{phyloseq-class} recommended.}

  \item{dist.fun}{(Optional). Default "jaccard". Any
  supported argument to the \code{method} parameter of the
  \code{\link{distance}} function is supported here. Some
  distance methods, like \code{"unifrac"}, may take a
  non-trivial amount of time to calculate, in which case
  you probably want to calculate the distance matrix
  separately, save, and then provide it as the argument to
  \code{dist.fun} instead. See below for alternatives).

  Alternatively, if you have already calculated the
  sample-wise distance object, the resulting
  \code{dist}-class object can be provided as
  \code{dist.fun} instead (see examples).

  A third alternative is to provide a function that takes a
  sample-by-species matrix (typical vegan orientation) and
  returns a sample-wise distance matrix.}

  \item{max.dist}{(Optional). Default \code{0.4}. The
  maximum ecological distance (as defined by
  \code{dist.fun}) allowed between two samples to still
  consider them ``connected'' by an edge in the graphical
  model.}

  \item{keep.isolates}{(Optional). Default \code{FALSE}.
  Logical. Whether to keep isolates (un-connected samples,
  not microbial isolates) in the graphical model that is
  returned. Default results in isolates being removed from
  the object.}

  \item{...}{(Optional). Additional parameters passed on to
  \code{\link{distance}}.}
}
\value{
  A \code{\link[igraph]{igraph}}-class object.
}
\description{
  A specialized function for creating graphical models of
  microbiome samples based on a user-defined ecological
  distance and threshold. The graph is ultimately built
  with tools from the \code{\link[igraph]{igraph}}-package.
}
\examples{
# # Example plots with Enterotype Dataset
data(enterotype)
ig <- make_sample_network(enterotype, max.dist=0.3)
plot_sample_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.3, label=NULL)
#
# ig <- make_sample_network(enterotype, max.dist=0.2)
# plot_sample_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.3, label=NULL)
#
# # Three methods of choosing/providing distance/distance-method
# Provide method name available to distance
ig <- make_sample_network(enterotype, max.dist=0.3, dist.fun="jaccard")
# Provide distance object, already computed
jaccdist <- distance(enterotype, "jaccard")
ih <- make_sample_network(enterotype, max.dist=0.3, dist.fun=jaccdist)
# Provide "custom" function.
ii <- make_sample_network(enterotype, max.dist=0.3, dist.fun=function(x){vegdist(x, "jaccard")})
# The have equal results:
all.equal(ig, ih)
all.equal(ig, ii)
#
# Try out making a trivial "network" of the 3-sample esophagus data,
# with weighted-UniFrac as distance
data(esophagus)
ij <- make_sample_network(esophagus, "unifrac", weighted=TRUE)
}
\seealso{
  \code{\link{plot_sample_network}}
}

