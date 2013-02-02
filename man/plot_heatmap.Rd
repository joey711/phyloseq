\name{plot_heatmap}
\alias{plot_heatmap}
\title{Create an ecologically-organized heatmap using ggplot2 graphics}
\usage{
  plot_heatmap(physeq, method="NMDS", distance="bray",
  sample.label=NULL, species.label=NULL, low="#000033",
  high="#66CCFF", na.value="black", trans=log_trans(4),
  max.label=250, title=NULL, ...)
}
\arguments{
  \item{physeq}{(Required). The data, in the form of an
  instance of the \code{\link{phyloseq-class}}. This should
  be what you get as a result from one of the
  \code{\link{import}} functions, or any of the processing
  downstream. No data components beyond the
  \code{\link{otu_table}} are strictly necessary, though
  they may be useful if you want to re-label the axis ticks
  according to some observable or taxonomic rank, for
  instance, or if you want to use a
  \code{\link{UniFrac}}-based distance (in which case your
  \code{physeq} data would need to have a tree included).}

  \item{method}{(Optional). The ordination method to use
  for organizing the heatmap. A great deal of the
  usefulness of a heatmap graphic depends upon the way in
  which the rows and columns are ordered.}

  \item{distance}{(Optional). A character string. The
  ecological distance method to use in the ordination. See
  \code{\link{distance}}.}

  \item{sample.label}{(Optional). A character string. The
  sample variable by which you want to re-label the sample
  (horizontal) axis.}

  \item{species.label}{(Optional). A character string. The
  taxonomic rank by which you want to re-label the
  species/OTU (vertical) axis.}

  \item{low}{(Optional). A character string. An R color.
  See \code{\link{colors}} for options support in R (there
  are lots). The color that represents the lowest non-zero
  value in the heatmap. Default is a dark blue color,
  \code{"#000033"}.}

  \item{high}{(Optional). A character string. An R color.
  See \code{\link{colors}} for options support in R (there
  are lots). The color that will represent the highest
  value in the heatmap. The default is \code{"#66CCFF"}.
  Zero-values are treated as \code{NA}, and set to
  \code{"black"}, to represent a background color.}

  \item{na.value}{(Optional). A character string. An R
  color. See \code{\link{colors}} for options support in R
  (there are lots). The color to represent what is
  essentially the background of the plot, the
  non-observations that occur as \code{NA} or \code{0}
  values in the abundance table. The default is
  \code{"black"}, which works well on computer-screen
  graphics devices, but may be a poor choice for printers,
  in which case you might want this value to be
  \code{"white"}, and reverse the values of \code{high} and
  \code{low}, above.}

  \item{trans}{(Optional). \code{"trans"}-class
  transformer-definition object. A numerical transformer to
  use in the continuous color scale. See
  \code{\link[scales]{trans_new}} for details. The default
  is \code{log_trans(4)}.}

  \item{max.label}{(Optional). Integer. Default is 250. The
  maximum number of labeles to fit on a given axis (either
  x or y). If number of taxa or samples exceeds this value,
  the corresponding axis will be stripped of any labels.

  This supercedes any arguments provided to
  \code{sample.label} or \code{species.label}. Make sure to
  increase this value if, for example, you want a special
  label for an axis that has 300 indices.}

  \item{title}{(Optional). Default \code{NULL}. Character
  string. The main title for the graphic.}

  \item{...}{(Optional). Additional parameters passed to
  \code{\link{ordinate}}.}
}
\value{
  A heatmap plot, in the form of a \code{\link{ggplot}}
  plot object, which can be further saved and modified.
}
\description{
  There are many useful examples of phyloseq heatmap
  graphics in the
  \href{http://joey711.github.com/phyloseq/plot_heatmap-examples}{phyloseq
  online tutorials}. In a 2010 article in BMC Genomics,
  Rajaram and Oono show describe an approach to creating a
  heatmap using ordination methods to organize the rows and
  columns instead of (hierarchical) cluster analysis. In
  many cases the ordination-based ordering does a much
  better job than h-clustering. An immediately useful
  example of their approach is provided in the NeatMap
  package for R. The NeatMap package can be used directly
  on the abundance table (\code{\link{otu_table-class}}) of
  phylogenetic-sequencing data, but the NMDS or PCA
  ordination options that it supports are not based on
  ecological distances. To fill this void, phyloseq
  provides the \code{plot_heatmap()} function as an
  ecology-oriented variant of the NeatMap approach to
  organizing a heatmap and build it using ggplot2 graphics
  tools. The \code{distance} and \code{method} arguments
  are the same as for the \code{\link{plot_ordination}}
  function, and support large number of distances and
  ordination methods, respectively, with a strong leaning
  toward ecology. This function also provides the options
  to re-label the OTU and sample axis-ticks with a
  taxonomic name and/or sample variable, respectively, in
  the hope that this might hasten your interpretation of
  the patterns (See the \code{sample.label} and
  \code{species.label} documentation, below). Note that
  this function makes no attempt to overlay hierarchical
  clustering trees on the axes, as hierarchical clustering
  is not used to organize the plot. Also note that each
  re-ordered axis repeats at the edge, and so apparent
  clusters at the far right/left or top/bottom of the
  heat-map may actually be the same. For now, the placement
  of this edge can be considered arbitrary, so beware of
  this artifact of this graphical representation. If you
  benefit from this phyloseq-specific implementation of the
  NeatMap approach, please cite both our packages/articles.
}
\details{
  This approach borrows heavily from the \code{heatmap1}
  function in the \code{NeatMap} package. Highly
  recommended, and we are grateful for their package and
  ideas, which we have adapted for our specific purposes
  here, but did not use an explicit dependency. At the time
  of the first version of this implementation, the NeatMap
  package depends on the rgl-package, which is not needed
  in phyloseq, at present. Although likely a transient
  issue, the rgl-package has some known installation issues
  that have further influenced to avoid making NeatMap a
  formal dependency (Although we love both NeatMap and
  rgl!).
}
\examples{
data("GlobalPatterns")
gpac <- subset_species(GlobalPatterns, Phylum=="Crenarchaeota")
# FYI, the base-R function uses a non-ecological ordering scheme,
# but does add potentially useful hclust dendrogram to the sides...
heatmap(otu_table(gpac))
plot_heatmap(gpac)
# example relabelling based on a sample variable and taxonomic rank.
plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family")
# Now repeat the plot, but change the color scheme in various ways.
# This may be worthwhile, depending on the graphics device or paper
# on which you want to discplay the heatmap.
plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family", low="#000033", high="#CCFF66")
plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family", low="#000033", high="#FF3300")
plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family", low="#000033", high="#66CCFF")
plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family", low="#66CCFF", high="#000033", na.value="white")
plot_heatmap(gpac, "NMDS", "bray", "SampleType", "Family", low="#FFFFCC", high="#000033", na.value="white")
# Now try the default color scheme, but using different ecological distances/ordinations
plot_heatmap(gpac, "NMDS", "jaccard")
plot_heatmap(gpac, "DCA")
plot_heatmap(gpac, "RDA")
plot_heatmap(gpac, "MDS", "bray")
plot_heatmap(gpac, "MDS", "unifrac")
plot_heatmap(gpac, "MDS", "unifrac", weighted=TRUE)
}
\references{
  Because this function relies so heavily in principle, and
  in code, on some of the functionality in NeatMap, please
  site their article if you use this function in your work.

  Rajaram, S., & Oono, Y. (2010). NeatMap--non-clustering
  heat map alternatives in R. BMC Bioinformatics, 11, 45.
}

