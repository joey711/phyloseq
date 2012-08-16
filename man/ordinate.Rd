\name{ordinate}
\alias{ordinate}
\title{Perform an ordination on phyloseq data}
\usage{
  ordinate(physeq, method="DCA", distance="unifrac", ...)
}
\arguments{
  \item{physeq}{(Required). Phylogenetic sequencing data
  (\code{\link{phyloseq-class}}). The data on which you
  want to perform the the ordination. In general, these
  methods will be based in some fashion on the abundance
  table ultimately stored as a contingency matrix
  (\code{\link{otu_table-class}}). If you're able to import
  data into \code{\link{phyloseq-class}} format, than you
  don't need to worry, as an \code{otu_table} is a required
  component of this class. In addition, some ordination
  methods require additional data, like a constraining
  variable or phylogenetic tree. If that is the case, the
  relevant data should be included in \code{physeq} prior
  to running. Integrating the data in this way also results
  in these different data components being checked for
  validity and completeness by the method.}

  \item{method}{(Optional). A character string. Default is
  \code{"DCA"}.

  Currently supported method options are: \code{c("DCA",
  "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")}

  \describe{ \item{DCA}{Performs detrended correspondence
  analysis using \code{\link{decorana}}}
  \item{CCA}{Performs correspondence analysis, or
  optionally, constrained correspondence analysis (a.k.a.
  canonical correspondence analysis), via
  \code{\link[vegan]{cca}} } \item{RDA}{Performs redundancy
  analysis, or optionally principal components analysis,
  via \code{\link[vegan]{rda}} } \item{DPCoA}{Performs
  Double Principle Coordinate Analysis using a (corrected,
  if necessary) phylogenetic/patristic distance between
  species. The calculation is performed by
  \code{\link{DPCoA}}(), which ultimately uses
  \code{\link[ade4]{dpcoa}} after making the appropriate
  accessions/corrections of the data. }
  \item{NMDS}{Performs Non-metric MultiDimenstional Scaling
  of a sample-wise ecological distance matrix onto a
  user-specified number of axes, \code{k}. By default,
  \code{k=2}, but this can be modified as a supplementary
  argument. This method is ultimately carried out by
  \code{\link{metaMDS}} after the appropriate accessions
  and distance calculations. Because \code{metaMDS}
  includes its own distance calculation wrappers to
  \code{\link[vegan]{vegdist}}, and these provide
  additional functionality in the form of species scores,
  \code{ordinate} will pass-on the \code{distance} argument
  to \code{metaMDS} if it is among the supported
  \code{vegdist} methods. However, all distance methods
  supported by \code{\link{distance}} are supported here,
  including \code{"unifrac"} (the default) and
  \code{"DPCoA"}. } \item{MDS/PCoA}{Performs principal
  coordinate analysis (also called principle coordinate
  decomposition, multidimensional scaling (MDS), or
  classical scaling) of a distance matrix (Gower 1966),
  including two correction methods for negative
  eigenvalues. See \code{\link[ape]{pcoa}} for further
  details. } }}

  \item{distance}{(Optional). A character string matching a
  \code{\link{distance}} method; or, alternatively, a
  pre-computed \code{\link{dist}}-class object. This
  argument is only utilized if a distance matrix is
  required by the ordination method specified by the
  \code{method} argument (above).

  Any supported \code{\link{distance}} methods are
  supported arguments to \code{distance} here. Try
  \code{distance("list")} for a explicitly supported
  distance method abbreviations. User-specified custom
  distance equations should also work, e.g.
  \code{"(A+B-2*J)/(A+B)"}. See \code{\link{distance}} for
  more details, examples.}

  \item{...}{(Optional). Additional arguments to supporting
  functions. For example, the additional argument
  \code{weighted=TRUE} would be passed on to
  \code{\link{UniFrac}} if \code{"unifrac"} were chosen as
  the \code{distance} option and \code{"MDS"} as the
  ordination \code{method} option. Alternatively, if
  \code{"DCA"} were chosen as the ordination \code{method}
  option, additional arguments would be passed on to the
  relevant ordination function, \code{\link{decorana}}, for
  example.}
}
\value{
  An ordination object. The specific class of the returned
  object depends upon the ordination method, as well as the
  function/package that is called internally to perform it.
  As a general rule, any of the ordination classes returned
  by this function will be recognized by downstream tools
  in the \code{phyloseq} package, for example the
  ordination plotting function,
  \code{\link{plot_ordination}}.
}
\description{
  This function wraps several commonly-used ordination
  methods. The type of ordination depends upon the argument
  to \code{method}. Try \code{ordinate("help")} or
  \code{ordinate("list")} for the currently supported
  method options.
}
\examples{
# # Take a subset of the GP dataset for quicker computation of examples
# data(GlobalPatterns)
# # Keep top 200 species
# topsp <- names(sort(taxa_sums(GlobalPatterns), TRUE)[1:200])
# GP    <- prune_taxa(topsp, GlobalPatterns)
# # Subset further to top 5 phyla
# top5ph <- sort(tapply(taxa_sums(GP), tax_table(GP)[, "Phylum"], sum), decreasing=TRUE)[1:5]
# GP     <- subset_taxa(GP, Phylum \%in\% names(top5ph))
# #
# # Examples performing ordination with NMDS. Default distance is unweighted UniFrac
# GP.NMDS.UF.ord   <- ordinate(GP, "NMDS")
# GP.NMDS.wUF.ord  <- ordinate(GP, "NMDS", "unifrac", weighted=TRUE)
# GP.NMDS.Bray.ord <- ordinate(GP, "NMDS", "bray")
# #
# # # An example plot with default, or manually-defined shapes
# (p <- plot_ordination(GP, GP.NMDS.Bray.ord, "biplot", color="SampleType", shape="Phylum"))
# # define manual shape scale:
# man.shapes <- 21:25
# names(man.shapes) <- c(get_taxa_unique(GP, "Phylum"))
# man.shapes <- c(samples=19, man.shapes)
# p + scale_shape_manual(value=man.shapes)
# #
# # An example of constrained ordination
# GP.cca <- ordinate(GP~SampleType, "CCA")
# #
# # Run-through "quick" plot examples of the other ordination options currently supported
# # Only showing "samples" in these examples, but "species" options supported for some methods
# plot_ordination(GP, ordinate(GP, "DCA"), "samples", color="SampleType")
# plot_ordination(GP, ordinate(GP, "CCA"), "samples", color="SampleType")
# plot_ordination(GP, ordinate(GP~SampleType, "CCA"), "samples", color="SampleType")
# plot_ordination(GP, ordinate(GP, "RDA"), "samples", color="SampleType")
# plot_ordination(GP, ordinate(GP~SampleType, "RDA"), "samples", color="SampleType")
# plot_ordination(GP, ordinate(GP, "DPCoA"), "samples", color="SampleType")
# plot_ordination(GP, ordinate(GP, "MDS"), "samples", color="SampleType")
# plot_ordination(GP, ordinate(GP, "PCoA"), "samples", color="SampleType")
# plot_ordination(GP, ordinate(GP, "NMDS"), "samples", color="SampleType")
# plot_ordination(GP, ordinate(GP, "NMDS", "w"), "samples", color="SampleType")
}
\seealso{
  Related component ordination functions described within
  phyloseq:

  \code{\link{DPCoA}}

  Described/provided by other packages:

  \code{\link{cca}}/\code{\link{rda}},
  \code{\link{decorana}}, \code{\link{metaMDS}},
  \code{\link{pcoa}}

  NMDS and MDS/PCoA both operate on distance matrices,
  typically based on some pairwise comparison of the
  microbiomes in an experiment/project. There are a number
  of common methods to use to calculate these pairwise
  distances, and the most convenient function (from a
  \code{phyloseq} point of view) for calculating these
  distance matrices is the

  \code{\link{distance}}

  function. It can be thought of as a distance /
  dissimilarity-index companion function for
  \code{ordinate}, and indeed the distance options provided
  to \code{ordinate} simply passed on to
  \code{\link{distance}}.

  A good quick summary of ordination is provided in the
  introductory vignette for vegan:

  \href{http://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf}{vegan
  introductory vignette}

  The following \code{R} task views are also useful for
  understanding the available tools in \code{R}:

  \href{http://cran.r-project.org/web/views/Environmetrics.html}{Analysis
  of Ecological and Environmental Data}

  \href{http://cran.r-project.org/web/views/Multivariate.html}{Multivariate
  Statistics}
}

