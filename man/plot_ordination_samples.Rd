\name{plot_ordination_samples}
\alias{plot_ordination_samples}
\title{Convenient rendering of ordination samples/sites using ggplot2.}
\usage{
  plot_ordination_samples(physeq, ordination, axes=c(1, 2),
  color=NULL, shape=NULL, label=NULL, title=NULL)
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

  \item{axes}{(Optional). A 2-element vector indicating the
  axes of the ordination that should be used for plotting.
  Can be \code{\link{character-class}} or
  \code{\link{integer-class}}, naming the index name or
  index of the desired axis for the horizontal and vertical
  axes, respectively, in that order. The default value,
  \code{c(1, 2)}, specifies the first two axes of the
  provided ordination.}

  \item{color}{(Optional). Default \code{NULL}. The sample
  variable to map to different colors. This can be a single
  character string of the variable name in
  \code{sampleData} (among the set returned by
  \code{sample.variables(physeq)} ); or a custom supplied
  vector with length equal to the number of samples in the
  dataset (nsamples(physeq)). The color scheme is chosen
  automatically by \code{link{ggplot}}, but it can be
  modified afterward with an additional layer using
  \code{\link[ggplot2]{scale_color_manual}}.}

  \item{shape}{(Optional). Default \code{NULL}. The sample
  variable to map to different shapes. Like \code{color},
  this can be a single character string of the variable
  name in \code{sampleData} (among the set returned by
  \code{sample.variables(physeq)} ); or a custom supplied
  vector with length equal to the number of samples in the
  dataset (nsamples(physeq)). The shape scale is chosen
  automatically by \code{link{ggplot}}, but it can be
  modified afterward with an additional layer using
  \code{\link[ggplot2]{scale_shape_manual}}.}

  \item{label}{(Optional). Default \code{NULL}. The sample
  variable to use for labelling each sample coordinate in
  the plot.}

  \item{title}{(Optional). Default \code{NULL}. Character
  string. The title to include over the plot.}
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
# data(GlobalPatterns)
# # Define a human-associated versus non-human binary variable:
# human.levels <- levels( getVariable(GlobalPatterns, "SampleType") ) \%in\%
#      c("Feces", "Mock", "Skin", "Tongue")
# human <- human.levels[getVariable(GlobalPatterns, "SampleType")]
# names(human) <- sample.names(GlobalPatterns)
# # Add the human variable to GlobalPatterns sample data
# sampleData(GlobalPatterns)$human <- human
# # Calculate the Bray-Curtis distance between each sample
# GP.dist <- vegdist(GlobalPatterns, "bray")
# # Perform principle coordinates analysis
# GP.dist.pcoa <- pcoa(GP.dist)
# (p1 <- plot_ordination_samples(GlobalPatterns, GP.dist.pcoa, shape=human, color="SampleType") )
}
\seealso{
  \code{\link{plot_ordination_biplot}}
}

