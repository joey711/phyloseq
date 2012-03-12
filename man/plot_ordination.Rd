\name{plot_ordination}
\alias{plot_ordination}
\title{General ordination plotter based on ggplot2.}
\usage{
  plot_ordination(physeq, ordination, type="samples",
  axes=c(1, 2), color=NULL, shape=NULL, label=NULL,
  title=NULL, justDF=FALSE)
}
\arguments{
  \item{physeq}{(Required). \code{\link{phyloseq-class}},
  or alternatively, an \code{\link{sampleData-class}}. The
  data about which you want to plot and annotate the
  ordination.}

  \item{ordination}{(Required). An ordination object. Many
  different classes of ordination are defined by \code{R}
  packages. The supported classes should be listed
  explicitly, but in the meantime, all ordination classes
  currently supported by the \code{\link[vegan]{scores}}
  function are supported here. There is no default, as the
  expectation is that the ordination will be performed and
  saved prior to calling this plot function.}

  \item{type}{(Optional). The plot type. Default is
  \code{"samples"}. The currently supported options are
  \code{c("samples", "sites", "species", "taxa", "biplot",
  "split")}. The option ``taxa'' is equivalent to
  ``species'' in this case, and similarly, ``samples'' is
  equivalent to ``sites''. The options \code{"sites"} and
  \code{"species"} result in a single-plot of just the
  sites/samples or species/taxa of the ordination,
  respectively. The \code{"biplot"} and \code{"split"}
  options result in a combined plot with both taxa and
  samples, either combined into one plot (``biplot'') or
  separated in two facet panels (``split''), respectively.}

  \item{axes}{(Optional). A 2-element vector indicating the
  axes of the ordination that should be used for plotting.
  Can be \code{\link{character-class}} or
  \code{\link{integer-class}}, naming the index name or
  index of the desired axis for the horizontal and vertical
  axes, respectively, in that order. The default value,
  \code{c(1, 2)}, specifies the first two axes of the
  provided ordination.}

  \item{color}{(Optional). Default \code{NULL}. Character
  string. The name of the variable to map to colors in the
  plot. This can be a sample variable (among the set
  returned by \code{sample.variables(physeq)} ) or
  taxonomic rank (among the set returned by
  \code{rank.names(physeq)}).

  Alternatively, if \code{type} indicates a single-plot
  (\code{"samples"} or \code{"species"}), then it is also
  possible to supply a custom vector with length equal to
  the relevant number of samples or species
  (\code{nsamples(physeq)} or \code{nspecies(physeq)}).

  Finally, The color scheme is chosen automatically by
  \code{link{ggplot}}, but it can be modified afterward
  with an additional layer using
  \code{\link[ggplot2]{scale_color_manual}}.}

  \item{shape}{(Optional). Default \code{NULL}. Character
  string. The name of the variable to map to different
  shapes on the plot. Similar to \code{color} option, but
  for the shape if points.

  The shape scale is chosen automatically by
  \code{link{ggplot}}, but it can be modified afterward
  with an additional layer using
  \code{\link[ggplot2]{scale_shape_manual}}.}

  \item{label}{(Optional). Default \code{NULL}. Character
  string. The name of the variable to map to text labels on
  the plot. Similar to \code{color} option, but for
  plotting text.}

  \item{title}{(Optional). Default \code{NULL}. Character
  string. The title to include over the plot.}

  \item{justDF}{(Optional). Default \code{FALSE}. Logical.
  Instead of returning a ggplot2-object, do you just want
  the relevant \code{data.frame} that was used to build the
  plot? This is a user-accessible option for obtaining the
  \code{data.frame}, in in principal to make a custom plot
  that isn't possible with the available options in this
  function. For contributing new functions (developers),
  the \code{\link{phyloseq-package}} provides/uses an
  internal function to build the key features of the
  \code{data.frame} prior to plot-build.}
}
\value{
  A \code{\link{ggplot}} plot object, graphically
  summarizing the ordination result for the specified axes.
}
\description{
  Convenience wrapper for plotting ordination results as a
  \code{ggplot2}-graphic, including additional annotation
  in the form of shading, shape, and/or labels of sample
  variables.
}
\examples{
##
# data(GlobalPatterns)
# # Define a human-associated versus non-human binary variable:
# human.levels <- levels( getVariable(GlobalPatterns, "SampleType") ) \%in\%
# c("Feces", "Mock", "Skin", "Tongue")
# human <- human.levels[getVariable(GlobalPatterns, "SampleType")]
# names(human) <- sample.names(GlobalPatterns)
# # Need to clean the zeros from GlobalPatterns:
# GP <- prune_species(speciesSums(GlobalPatterns)>0, GlobalPatterns)
# # Get the names of the most-abundant
# top.TaxaGroup <- sort(
# tapply(speciesSums(GP), taxTab(GP)[, "Phylum"], sum, na.rm = TRUE),
# decreasing = TRUE)
# top.TaxaGroup <- top.TaxaGroup[top.TaxaGroup > 1*10^6]
# # Now prune further, to just the most-abundant phyla
# GP <- subset_species(GP, Phylum \%in\% names(top.TaxaGroup))
# # DPCoA - because scores.dpcoa() is internal to phyloseq, can only test
# # this chunk when new plot_ordination() is installed
# topsp <- names(sort(speciesSums(GP), TRUE)[1:200])
# GP1   <- prune_species(topsp, GP)
# GP.dpcoa <- DPCoA(GP1)
# plot_ordination(GP1, GP.dpcoa, type="taxa", color="Phylum")
# plot_ordination(GP1, GP.dpcoa, type="samples", color="SampleType") + geom_line() + geom_point(size=5)
# plot_ordination(GP1, GP.dpcoa, type="samples", color="SampleType", shape=human) +
# geom_line() + geom_point(size=5)
# plot_ordination(GP1, GP.dpcoa, type="species", color="Phylum") + geom_line() + geom_point(size=5)
# plot_ordination(GP1, GP.dpcoa, type="biplot", shape="Phylum", label="SampleType")
# plot_ordination(GP1, GP.dpcoa, type="biplot", shape="Phylum")
# plot_ordination(GP1, GP.dpcoa, type="biplot", color="Phylum")
# plot_ordination(GP1, GP.dpcoa, type="biplot", label="Phylum")
# plot_ordination(GP1, GP.dpcoa, type="split", color="Phylum", label="SampleType")
# plot_ordination(GP1, GP.dpcoa, type="split", color="SampleType", shape="Phylum", label="SampleType")
}
\seealso{
  \code{\link{plot_phyloseq}}
}

