\name{subset_ord_plot}
\alias{subset_ord_plot}
\title{Subset points from an ordination-derived ggplot}
\usage{
  subset_ord_plot(p, threshold=0.05, method="farthest")
}
\arguments{
  \item{p}{(Required).  A \code{\link{ggplot}} object
  created by \code{\link{plot_ordination}}. It contains the
  complete data that you want to subset.}

  \item{threshold}{(Optional). A numeric scalar. Default is
  \code{0.05}. This value determines a coordinate threshold
  or population threshold, depending on the value of the
  \code{method} argument, ultimately determining which
  points are included in returned \code{data.frame}.}

  \item{method}{(Optional). A character string. One of
  \code{c("farthest", "radial", "square")}. Default is
  \code{"farthest"}. This determines how threshold will be
  interpreted.

  \describe{

  \item{farthest}{ Unlike the other two options, this
  option implies removing a certain fraction or number of
  points from the plot, depending on the value of
  \code{threshold}. If \code{threshold} is greater than or
  equal to \code{1}, then all but \code{threshold} number
  of points farthest from the origin are removed.
  Otherwise, if \code{threshold} is less than \code{1}, all
  but \code{threshold} fraction of points farthests from
  origin are retained. }

  \item{radial}{ Keep only those points that are beyond
  \code{threshold} radial distance from the origin. Has the
  effect of removing a circle of points from the plot,
  centered at the origin. }

  \item{square}{ Keep only those points with at least one
  coordinate greater than \code{threshold}. Has the effect
  of removing a ``square'' of points from the plot,
  centered at the origin. }

  }}
}
\value{
  A \code{\link{data.frame}} suitable for creating a
  \code{\link{ggplot}} plot object, graphically summarizing
  the ordination result according to previously-specified
  parameters.
}
\description{
  Easily retrieve a plot-derived \code{data.frame} with a
  subset of points according to a threshold and method. The
  meaning of the threshold depends upon the method. See
  argument description below.
}
\examples{
##
# data(GlobalPatterns)
# # Need to clean the zeros from GlobalPatterns:
# GP <- prune_species(speciesSums(GlobalPatterns)>0, GlobalPatterns)
# sampleData(GP)$human <- factor(human)
# # Get the names of the most-abundant phyla
# top.TaxaGroup <- sort(
#   	tapply(speciesSums(GP), taxTab(GP)[, "Phylum"], sum, na.rm = TRUE),
#   decreasing = TRUE)
# top.TaxaGroup <- top.TaxaGroup[top.TaxaGroup > 1*10^6]
# # Prune to just the most-abundant phyla
# GP <- subset_species(GP, Phylum \%in\% names(top.TaxaGroup))
# # Perform a correspondence analysis
# gpca <- ordinate(GP, "CCA")
# # # Make species topo with a subset of points layered
# # First, make a basic plot of just the species
# p1 <- plot_ordination(GP, gpca, "species", color="Phylum")
# # Re-draw this as topo without points, and facet
# p1 <- ggplot(p1$data, p1$mapping) + geom_density2d() + facet_wrap(~Phylum)
# # Add a layer of a subset of species-points that are furthest from origin.
# p53 <- p1 + geom_point(data=subset_ord_plot(p1, 1.0, "square"), size=1)
# print(p53)
}
\seealso{
  \code{\link{plot_ordination}}
}

