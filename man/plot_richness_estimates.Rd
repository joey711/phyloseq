\name{plot_richness_estimates}
\alias{plot_richness_estimates}
\title{Plot richness estimates, flexibly with ggplot2}
\usage{
  plot_richness_estimates(physeq, x, color=NULL,
  shape=NULL)
}
\arguments{
  \item{physeq}{(Required). \code{\link{phyloseq-class}},
  or alternatively, an \code{\link{otuTable-class}}. The
  data about which you want to estimate the richness.}

  \item{x}{(Optional). A variable to map to the horizontal
  axis. The vertical axis will be mapped to richness
  estimates and have units of total species. This parameter
  (\code{x}) can be either a character string indicating a
  variable in \code{sampleData} (among the set returned by
  \code{sample.variables(physeq)} ); or a custom supplied
  vector with length equal to the number of samples in the
  dataset (nsamples(physeq)).

  The default value is \code{"sample.names"}, which will
  map each sample's name to a separate horizontal position
  in the plot.}

  \item{color}{(Optional). Default \code{NULL}. The sample
  variable to map to different colors. Like \code{x}, this
  can be a single character string of the variable name in
  \code{sampleData} (among the set returned by
  \code{sample.variables(physeq)} ); or a custom supplied
  vector with length equal to the number of samples in the
  dataset (nsamples(physeq)). The color scheme is chosen
  automatically by \code{link{ggplot}}, but it can be
  modified afterward with an additional layer using
  \code{\link[ggplot2]{scale_color_manual}}.}

  \item{shape}{(Optional). Default \code{NULL}. The sample
  variable to map to different shapes. Like \code{x} and
  \code{color}, this can be a single character string of
  the variable name in \code{sampleData} (among the set
  returned by \code{sample.variables(physeq)} ); or a
  custom supplied vector with length equal to the number of
  samples in the dataset (nsamples(physeq)). The shape
  scale is chosen automatically by \code{link{ggplot}}, but
  it can be modified afterward with an additional layer
  using \code{\link[ggplot2]{scale_shape_manual}}.}
}
\value{
  A \code{\link{ggplot}} plot object summarizing the
  richness estimates, and their standard error.
}
\description{
  Performs a number of standard richness estimates using
  the \code{\link{estimate_richness}} function, and returns
  a \code{ggplot} plotting object. This plot shows the
  individual richness estimates for each sample, as well as
  the observed richness. You must use untrimmed datasets
  for meaningful results, as these estimates (and even the
  ``observed'' richness) are highly dependent on the number
  of singletons. You can always trim the data later on if
  needed, just not before using this function.
}
\details{
  NOTE: Because this plotting function incorporates the
  output from \code{\link{estimate_richness}}, the variable
  names of that output should not be used as \code{x} or
  \code{color} (even if it works, the resulting plot might
  be kindof strange, and not the intended behavior of this
  function). The following are the names you will want to
  avoid using in \code{x} or \code{color}:

  \code{c("S.obs", "S.chao1", "se.chao1", "S.ACE",
  "se.ACE", "shannon", "simpson")}
}
\examples{
# data(GlobalPatterns)
# plot_richness_estimates(GlobalPatterns, "SampleType")
# plot_richness_estimates(GlobalPatterns, "SampleType", "SampleType")
#
# # Define a human-associated versus non-human categorical variable:
# GP <- GlobalPatterns
# human.levels <- levels( getVariable(GP, "SampleType") ) \%in\%
# c("Feces", "Mock", "Skin", "Tongue")
# human <- human.levels[getVariable(GP, "SampleType")]
# names(human) <- sample.names(GP)
# # Replace current SD with new one that includes human variable:
# sampleData(GP) <- sampleData(data.frame(sampleData(GP), human))
#
# # Can use new "human" variable within GP as a discrete variable in the plot
# plot_richness_estimates(GP, "human", "SampleType")
# plot_richness_estimates(GP, "SampleType", "human")
#
# # Can also provide custom factor directly:
# plot_richness_estimates(GP, "SampleType", human)
# plot_richness_estimates(GP, human, "SampleType")
#
# # Not run: Should cause an error:
# plot_richness_estimates(GP, "value", "value")
# #
}
\seealso{
  \code{\link{estimate_richness}},
  \code{\link[vegan]{estimateR}},
  \code{\link[vegan]{diversity}}
}

