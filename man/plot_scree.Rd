\name{plot_scree}
\alias{plot_scree}
\title{General ordination eigenvalue plotter using ggplot2.}
\usage{
  plot_scree(ordination, title = NULL)
}
\arguments{
  \item{ordination}{(Required). An ordination object. Many
  different classes of ordination are defined by \code{R}
  packages. Ordination classes currently supported/created
  by the \code{\link{ordinate}} function are supported
  here. There is no default, as the expectation is that the
  ordination will be performed and saved prior to calling
  this plot function.}

  \item{title}{(Optional). Default \code{NULL}. Character
  string. The main title for the graphic.}
}
\value{
  A \code{\link{ggplot}} plot object, graphically
  summarizing the ordination result for the specified axes.
}
\description{
  Convenience wrapper for plotting ordination eigenvalues
  (if available) using a \code{ggplot2}-graphic.
}
\examples{
# First load and trim a dataset
data(GlobalPatterns)
GP = prune_taxa(taxa_sums(GlobalPatterns)>0, GlobalPatterns)
# Define a human-associated versus non-human categorical variable, and add new human variable to sample data:
sample_data(GP)$human = factor( get_variable(GP, "SampleType") \%in\% c("Feces", "Mock", "Skin", "Tongue") )
# # filtering
# Remove taxa not seen more than 3 times in at least 20\% of the samples
gp  = filter_taxa(GP, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
# Standardize abundances to the median sequencing depth
gpr = transform_sample_counts(gp, function(x, total=median(sample_sums(gp))) round(total * (x / sum(x))) )
# Let's use Coefficient of Variation for filtering, arbitrary cutoff of 3.0
gprf = filter_taxa(gpr, function(x) sd(x)/mean(x) > 3L, TRUE)
# For a somewhat readable number of taxa on display, let's subset to just Bacteroidetes for some plots
gprfb = subset_taxa(gprf, Phylum=="Bacteroidetes")
# Test plots (preforms ordination in-line, then makes scree plot)
plot_scree(ordinate(gprfb, "DPCoA", "bray"))
plot_scree(ordinate(gprfb, "PCoA", "bray"))
plot_scree(ordinate(gprfb, "NMDS", "bray")) # Empty return with message
plot_scree(ordinate(gprfb ~ SampleType, "CCA"))
plot_scree(ordinate(gprfb ~ SampleType, "RDA"))
plot_scree(ordinate(gprfb, "DCA"))
plot_ordination(gprfb, ordinate(gprfb, "DCA"), type="scree")
}
\seealso{
  \code{\link{plot_ordination}}

  \code{\link{ordinate}}

  \code{\link{distance}}

  The examples on the phyloseq wiki page for
  \code{plot_ordination} show many more examples:

  \url{https://github.com/joey711/phyloseq/wiki/plot_ordination}
}

