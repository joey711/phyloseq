\name{plot_tree}
\alias{plot_tree}
\title{Plot a phylogenetic tree with optional annotations}
\usage{
  plot_tree(physeq, method="sampledodge", color=NULL,
  shape=NULL, size=NULL, min.abundance=Inf,
  label.tips=NULL, text.size=NULL, sizebase=5, base.spacing
  = 0.02)
}
\arguments{
  \item{physeq}{(Required). The data about which you want
  to plot and annotate a phylogenetic tree, in the form of
  a single instance of the \code{\link{phyloseq-class}},
  containing at minimum a phylogenetic tree component (try
  \code{\link{tre}}). One of the major advantages of this
  function over basic tree-plotting utilities in the
  \code{\link{ape}}-package is the ability to easily
  annotate the tree with sample variables and taxonomic
  information. For these uses, the \code{physeq} argument
  should also have a \code{\link{sampleData}} and/or
  \code{\link{taxTab}} component(s).}

  \item{method}{(Optional). Character string. Default
  \code{"sampledodge"}. The name of the annotation method
  to use. This will be expanded in future versions.
  Currently only \code{"sampledodge"} and \code{"treeonly"}
  are supported. The \code{"sampledodge"} option results in
  points drawn next to leaves if individuals from that taxa
  were observed, and a separate point is drawn for each
  sample.}

  \item{color}{(Optional). Character string. Default
  \code{NULL}. The name of the variable in \code{physeq} to
  map to point color.}

  \item{shape}{(Optional). Character string. Default
  \code{NULL}. The name of the variable in \code{physeq} to
  map to point shape.}

  \item{size}{(Optional). Character string. Default
  \code{NULL}. The name of the variable in \code{physeq} to
  map to point size. A special argument \code{"abundance"}
  is reserved here and scales point size using abundance in
  each sample on a log scale.}

  \item{min.abundance}{(Optional). Numeric. The minimum
  number of individuals required to label a point with the
  precise number. Default is \code{Inf}, meaning that no
  points will have their abundance labeled. If a vector,
  only the first element is used.}

  \item{label.tips}{(Optional). Character string. Default
  is \code{NULL}, indicating that no tip labels will be
  printed. If \code{"species.names"}, then the name of the
  taxa will be added to the tree; either next to the
  leaves, or next to the set of points that label the
  leaves. Alternatively, if this is one of the rank names
  (from \code{rank.names(physeq)}), then the identity (if
  any) for that particular taxonomic rank is printed
  instead.}

  \item{text.size}{(Optional). Numeric. Should be positive.
  The size parameter used to control the text size of taxa
  labels. Default is \code{NULL}. If left \code{NULL}, this
  function will automatically calculate a (hopefully)
  optimal text size given the vertical constraints posed by
  the tree itself. This argument is included in case the
  automatically-calculated size is wrong, and you want to
  change it. Note that this parameter is only meaningful if
  \code{label.tips} is not \code{NULL}.}

  \item{sizebase}{(Optional). Numeric. Should be positive.
  The base of the logarithm used to scale point sizes to
  graphically represent abundance of species in a given
  sample. Default is 5.}

  \item{base.spacing}{(Optional). Numeric. Default is
  \code{0.02}. Should be positive. This defines the
  base-spacing between points at each tip/leaf in the the
  tree. The larger this value, the larger the spacing
  between points. This is useful if you have problems with
  overlapping large points and/or text indicating
  abundance, for example. Similarly, if you don't have this
  problem and want tighter point-spacing, you can shrink
  this value.}
}
\value{
  A \code{\link{ggplot}}2 plot.
}
\description{
  This function is intended to facilitate easy graphical
  investigation of the phylogenetic tree, as well as sample
  data. Note that for phylogenetic sequencing of samples
  with large richness, some of the options in this function
  will be prohibitively slow to render, or too dense to be
  interpretable. A rough ``rule of thumb'' is to use
  subsets of data with not many more than 200 taxa per
  plot, sometimes less depending on the complexity of the
  additional annotations being mapped to the tree. It is
  usually possible to create an unreadable, uninterpretable
  tree with modern datasets. However, the goal should be
  toward parameter settings and data subsets that convey
  (honestly, accurately) some biologically relevant feature
  of the data. One of the goals of the
  \code{\link{phyloseq-package}} is to make the
  determination of these features/settings as easy as
  possible.
}
\details{
  This function received a development contribution from
  the work of Gregory Jordan for the \code{ggphylo}
  package, which provides tools for rendering a
  phylogenetic tree in \code{\link{ggplot2}} graphics. That
  package is not currently available from CRAN or
  Bioconductor, but is available in development (roughly
  ``beta'') form from GitHub. Furthermore, although
  \code{ggphylo} awesomely supports radial and unrooted
  trees, \code{plot_tree} currently only supports
  ``standard'' square-horizontal trees. Additional support
  for these types of features (like radial trees) is
  planned. Send us development feedback if this is a
  feature you really want to have soon.
}
\examples{
# # Using plot_tree() with the esophagus dataset.
data(esophagus)
plot_tree(esophagus)
plot_tree(esophagus, color="samples")
plot_tree(esophagus, size="abundance")
plot_tree(esophagus, size="abundance", color="samples")
plot_tree(esophagus, size="abundance", color="samples", base.spacing=0.03)
# # Using plot_tree with the Global Patterns dataset
data("GlobalPatterns")
# Subset Global Patterns dataset to just the observed Archaea
gpa <- subset_species(GlobalPatterns, Kingdom=="Archaea")
# The number of different Archaeal species from this dataset is small enough ...
nspecies(gpa)
# ... that it is reasonable to consider displaying the phylogenetic tree directly.
# (probably not true of the total dataset)
nspecies(GlobalPatterns)
# Some patterns are immediately discernable with minimal parameter choices:
# plot_tree(gpa, color="SampleType")
# plot_tree(gpa, color="Phylum")
# plot_tree(gpa, color="SampleType", shape="Phylum")
# plot_tree(gpa, color="Phylum", label.tips="Genus")
# # However, the text-label size scales with number of species, and with common
# # graphics-divice sizes/resolutions, these ~200 taxa still make for a
# # somewhat crowded graphic.
# #
# # Let's instead subset further ot just the Crenarchaeota
gpac <- subset_species(gpa, Phylum=="Crenarchaeota")
plot_tree(gpac, color="SampleType", shape="Genus")
# plot_tree(gpac, color="SampleType", label.tips="Genus")
# # Let's add some abundance information.
# # Notice that the default spacing gets a little crowded when we map
# # species-abundance to point-size:
# plot_tree(gpac, color="SampleType", shape="Genus", size="abundance")
# # So let's spread it out a little bit with the base.spacing parameter.
# plot_tree(gpac, color="SampleType", shape="Genus", size="abundance", base.spacing=0.05)
}
\author{
  Paul McMurdie, relying on supporting code from Gregory
  Jordan \email{gjuggler@gmail.com}
}
\seealso{
  \code{\link{plot.phylo}}

  This function is a special use-case that relies on code
  borrowed directly, with permission, from the
  not-yet-officially-released package, \code{ggphylo},
  currently only available from GitHub at:
  \url{https://github.com/gjuggler/ggphylo}
}

