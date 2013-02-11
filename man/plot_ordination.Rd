\name{plot_ordination}
\alias{plot_ordination}
\title{General ordination plotter based on ggplot2.}
\usage{
  plot_ordination(physeq, ordination, type="samples",
  axes=c(1, 2), color=NULL, shape=NULL, label=NULL,
  title=NULL, justDF=FALSE)
}
\arguments{
  \item{physeq}{(Required). \code{\link{phyloseq-class}}.
  The data about which you want to plot and annotate the
  ordination.}

  \item{ordination}{(Required). An ordination object. Many
  different classes of ordination are defined by \code{R}
  packages. Ordination classes currently supported/created
  by the \code{\link{ordinate}} function are supported
  here. There is no default, as the expectation is that the
  ordination will be performed and saved prior to calling
  this plot function.}

  \item{type}{(Optional). The plot type. Default is
  \code{"samples"}. The currently supported options are
  \code{c("samples", "sites", "species", "taxa", "biplot",
  "split", "scree")}. The option ``taxa'' is equivalent to
  ``species'' in this case, and similarly, ``samples'' is
  equivalent to ``sites''. The options \code{"sites"} and
  \code{"species"} result in a single-plot of just the
  sites/samples or species/taxa of the ordination,
  respectively. The \code{"biplot"} and \code{"split"}
  options result in a combined plot with both taxa and
  samples, either combined into one plot (``biplot'') or
  separated in two facet panels (``split''), respectively.
  The \code{"scree"} option results in a call to
  \code{\link{plot_scree}}, which produces an ordered bar
  plot of the normalized eigenvalues associated with each
  ordination axis.}

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
  returned by \code{sample_variables(physeq)} ) or
  taxonomic rank (among the set returned by
  \code{rank_names(physeq)}).

  Alternatively, if \code{type} indicates a single-plot
  (\code{"samples"} or \code{"species"}), then it is also
  possible to supply a custom vector with length equal to
  the relevant number of samples or species
  (\code{nsamples(physeq)} or \code{ntaxa(physeq)}).

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
  string. The main title for the graphic.}

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
  There are many useful examples of phyloseq ordination
  graphics in the
  \href{http://joey711.github.com/phyloseq/plot_ordination-examples}{phyloseq
  online tutorials}. Convenience wrapper for plotting
  ordination results as a \code{ggplot2}-graphic, including
  additional annotation in the form of shading, shape,
  and/or labels of sample variables.
}
\examples{
data(GlobalPatterns)
# Need to clean the zeros from GlobalPatterns:
GP <- prune_taxa(taxa_sums(GlobalPatterns)>0, GlobalPatterns)
# Define a human-associated versus non-human binary variable:
sample_data(GP)$human <- get_variable(GP, "SampleType") \%in\% c("Feces", "Mock", "Skin", "Tongue")
# Get the names of the most-abundant
top.TaxaGroup <- sort(
tapply(taxa_sums(GP), tax_table(GP)[, "Phylum"], sum, na.rm = TRUE),
decreasing = TRUE)
top.TaxaGroup <- top.TaxaGroup[top.TaxaGroup > 1*10^6]
# Now prune further, to just the most-abundant phyla
GP <- subset_taxa(GP, Phylum \%in\% names(top.TaxaGroup))
topsp <- names(sort(taxa_sums(GP), TRUE)[1:200])
GP1   <- prune_taxa(topsp, GP)
GP.dpcoa <- ordinate(GP1, "DPCoA")
plot_ordination(GP1, GP.dpcoa, type="taxa", color="Phylum")
# Customize with ggplot2 layers added directly to output
library("ggplot2")
plot_ordination(GP1, GP.dpcoa, type="samples", color="SampleType") + geom_line() + geom_point(size=5)
p <- plot_ordination(GP1, GP.dpcoa, type="samples", color="SampleType", shape="human")
print(p)
# library("ggplot2")
# p + geom_line() + geom_point(size=5)
# plot_ordination(GP1, GP.dpcoa, type="taxa", color="Phylum") + geom_line() + geom_point(size=5)
plot_ordination(GP1, GP.dpcoa, type="biplot", shape="Phylum", label="SampleType")
plot_ordination(GP1, GP.dpcoa, type="biplot", shape="Phylum")
plot_ordination(GP1, GP.dpcoa, type="biplot", color="Phylum")
plot_ordination(GP1, GP.dpcoa, type="biplot", label="Phylum")
plot_ordination(GP1, GP.dpcoa, type="split", color="Phylum", label="SampleType")
plot_ordination(GP1, GP.dpcoa, type="split", color="SampleType", shape="Phylum", label="SampleType")
plot_ordination(GP1, GP.dpcoa, type="scree")
}
\seealso{
  The examples on the phyloseq wiki page for
  \code{plot_ordination} show many more examples:

  \url{https://github.com/joey711/phyloseq/wiki/plot_ordination}

  Also see the general wrapping function:

  \code{\link{plot_phyloseq}}
}

