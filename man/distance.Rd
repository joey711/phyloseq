\name{distance}
\alias{distance}
\title{General distance / dissimilarity index calculator}
\usage{
  distance(physeq, method="unifrac", type="samples", ...)
}
\arguments{
  \item{physeq}{(Required).  A \code{\link{phyloseq-class}}
  or an \code{\link{otuTable-class}} object. The latter is
  only appropriate for methods that do not require any
  additional data (one-table). For example, the ``unifrac''
  option (\code{\link{UniFrac}}) requires
  \code{\link{phyloseq-class}} that contains both an
  \code{otuTable} and a phylogenetic tree (\code{phylo}).}

  \item{method}{(Optional). A character string. Provide one
  of the 43 currently supported options. Default is
  \code{"unifrac"}. Alternatively, you can provide a
  character string that defines a custom distance method,
  if it has the form described in \code{\link{designdist}}.
  Depending on the \code{method} argument, this function
  wraps one of \code{\link{UniFrac}}, \code{\link{DPCoA}},
  \code{\link{vegdist}}, \code{\link{betadiver}},
  \code{\link{designdist}}, or \code{\link{dist}}.}

  \item{type}{(Optional). A character string. The type of
  pairwise comparisons being calculated: sample-wise or
  species-wise. The default is \code{c("samples")}.}

  \item{...}{Additional arguments passed on to the
  appropriate distance function, determined by the
  \code{method} argument.}
}
\value{
  An object of class ``\code{\link{dist}}'' suitable for
  certain ordination methods and other distance-based
  analyses.
}
\description{
  Takes a \code{\link{phyloseq-class}} object and method
  option, and returns a \code{\link{dist}}ance object
  suitable for certain ordination methods and other
  distance-based analyses. The are currently 43 explicitly
  supported method options, as well as user-provided
  arbitrary methods via an interface to
  \code{\link{designdist}}. Currrently only sample-wise
  distances are supported (the \code{type} argument), but
  eventually species-wise (OTU-wise) distances will be
  supported as well.
}
\examples{
data(esophagus)
distance(esophagus) # Unweighted UniFrac
distance(esophagus, weighted=TRUE) # weighted UniFrac
distance(esophagus, "jaccard") # vegdist jaccard
distance(esophagus, "gower") # vegdist option "gower"
distance(esophagus, "g") # designdist method option "g"
distance(esophagus, "minkowski") # invokes a method from the base dist() function.
distance(esophagus, "(A+B-2*J)/(A+B)") # designdist custom distance
distance("help")
distance("list")
help("distance")
}
\seealso{
  \code{\link{plot_ordination}}
}

