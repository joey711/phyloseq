\name{make_sample_network}
\alias{make_sample_network}
\title{Make sample-wise microbiome network (igraph)}
\usage{
  make_sample_network(physeq, keep.isolates=FALSE,
  threshold = 0, max.dist = 0.4, dist.fun = function(x)
  vegdist(x, "jaccard"), presence.absence=FALSE)
}
\arguments{
  \item{physeq}{(Required). Default \code{NULL}. A
  \code{\link{phyloseq-class}} object on which \code{g} is
  based.}

  \item{keep.isolates}{(Optional). Default \code{FALSE}.
  Logical. Whether to keep isolates (un-connected samples,
  not microbial isolates) in the graphical model that is
  returned. Default results in isolates being removed from
  the object.}

  \item{threshold}{(Optional). Default \code{1}. The
  minimum number of individuals of each taxa across all
  samples. Taxa below this value are ignored during network
  construction. Note that this is extremely simple
  filtering. More advanced filtering can be performed prior
  to this function using the \code{\link{filterfunSample}}
  or \code{\link{filter_taxa}} functions.}

  \item{max.dist}{(Optional). Default \code{0.4}. The
  maximum ecological distance (as defined by
  \code{dist.fun}) allowed between two samples to still
  consider them ``connected'' by an edge in the graphical
  model.}

  \item{dist.fun}{(Optional). Default \code{function(x)
  vegdist(x, "jaccard")}. A function that takes a
  \code{\link{phyloseq-class}} object (or
  \code{\link{otuTable-class}}) and returns a sample-wise
  distance matrix. Currently, any method in
  \code{\link{vegdist}} is supported, as well as
  \code{\link{UniFrac}} (although this may take some time
  to run you probably want to run this separately, save,
  and then provide the distance-object directly. See below
  for alternatives). Most ecological distance functions
  that take a matrix in \code{vegan-package} format
  (samples are rows) are supported.

  Alternatively, if you have already calculated the
  sample-wise distance object, this can be provided as
  \code{dist.fun} instead (see examples). A third
  alternative is to provide a character string indicating
  the desired method of distance calculation. For the
  moment, this is limited to the distance methods available
  in \code{\link{vegdist}}.}

  \item{presence.absence}{(Optional). Default \code{FALSE}.
  Whether the abundances values should be transformed to
  binary presence-absence values prior to distance
  calculation and network construction.}
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

ig <- make_sample_network(enterotype, FALSE, max.dist=0.3)
plot_sample_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.3, label=NULL)

# ig <- make_sample_network(enterotype, FALSE, max.dist=0.2)
# plot_sample_network(ig, enterotype, color="SeqTech", shape="Enterotype", line_weight=0.3, label=NULL)

# # Three methods of choosing/providing distance/distance-method
# Provide method name available to vegdist
ig <- make_sample_network(enterotype, FALSE, max.dist=0.3, dist.fun="jaccard")
# Provide distance object, already computed
ih <- make_sample_network(enterotype, FALSE, max.dist=0.3, dist.fun=vegdist(enterotype, "jaccard"))
# Provide explicit function. Can be custom function.
ii <- make_sample_network(enterotype, FALSE, max.dist=0.3, dist.fun=function(x){vegdist(x, "jaccard")})
# The have equal results:
all.equal(ig, ih)
all.equal(ig, ii)
#
}
\seealso{
  \code{\link{plot_sample_network}}
}

