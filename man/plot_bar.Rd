\name{plot_bar}
\alias{plot_bar}
\title{A flexible, informative barplot phyloseq data}
\usage{
  plot_bar(physeq, x="Sample", y="Abundance", fill=NULL,
  title=NULL, facet_grid=NULL)
}
\arguments{
  \item{physeq}{(Required). An
  \code{\link{otu_table-class}} or
  \code{\link{phyloseq-class}}.}

  \item{x}{(Optional). Optional, but recommended,
  especially if your data is comprised of many samples. A
  character string. The variable in the melted-data that
  should be mapped to the x-axis. See \code{\link{psmelt}},
  \code{\link{melt}}, and \code{\link{ggplot}} for more
  details.}

  \item{y}{(Optional). A character string. The variable in
  the melted-data that should be mapped to the y-axis.
  Typically this will be \code{"Abundance"}, in order to
  quantitatively display the abundance values for each
  OTU/group. However, alternative variables could be used
  instead, producing a very different, though possibly
  still informative, plot. See \code{\link{psmelt}},
  \code{\link{melt}}, and \code{\link{ggplot}} for more
  details.}

  \item{fill}{(Optional). A character string. Indicates
  which sample variable should be used to map to the fill
  color of the bars. The default is \code{NULL}, resulting
  in a gray fill for all bar segments.}

  \item{facet_grid}{(Optional). A formula object. It should
  describe the faceting you want in exactly the same way as
  for \code{\link[ggplot2]{facet_grid}}, and is ulitmately
  provided to \code{\link{ggplot}}2 graphics. The default
  is: \code{NULL}, resulting in no faceting.}

  \item{title}{(Optional). Default \code{NULL}. Character
  string. The main title for the graphic.}
}
\value{
  A \code{\link[ggplot2]{ggplot}}2 graphic object --
  rendered in the graphical device as the default
  \code{\link[base]{print}}/\code{\link[methods]{show}}
  method.
}
\description{
  There are many useful examples of phyloseq barplot
  graphics in the
  \href{http://joey711.github.com/phyloseq/plot_bar-examples}{phyloseq
  online tutorials}. This function wraps \code{ggplot2}
  plotting, and returns a \code{ggplot2} graphic object
  that can be saved or further modified with additional
  layers, options, etc. The main purpose of this function
  is to quickly and easily create informative summary
  graphics of the differences in taxa abundance between
  samples in an experiment.
}
\examples{
data("GlobalPatterns")
gp.ch = subset_species(GlobalPatterns, Phylum == "Chlamydiae")
plot_bar(gp.ch)
plot_bar(gp.ch, fill="Genus")
plot_bar(gp.ch, x="SampleType", fill="Genus")
plot_bar(gp.ch, "SampleType", fill="Genus", facet_grid=~Family)
# Need to load ggplot2 for custom faceting, and other custom ggplot2 goodies.
library("ggplot2")
plot_bar(gp.ch, fill="Genus") + facet_wrap(~SampleType)
plot_bar(gp.ch, fill="Genus") + facet_grid(SampleType ~ Family)
plot_bar(gp.ch, "SampleType", fill="Genus", facet_grid=SampleType~Family)
# A more complicated addition. Two lines. Second adds abundance points.
p = plot_bar(gp.ch, "SampleType", fill="Genus", facet_grid=~Genus)
p + geom_point(aes(x=SampleType, y=Abundance), color="black", position="jitter", size=1.5)
# Enterotype Example
data("enterotype")
TopNOTUs <- names(sort(taxa_sums(enterotype), TRUE)[1:10])
ent10   <- prune_species(TopNOTUs, enterotype)
plot_bar(ent10, "SeqTech", fill="Enterotype", facet_grid=~Genus)
# The previous was probably more informative, but here is the same
# information presented with a different organization.
plot_bar(ent10, "SeqTech", fill="Genus", facet_grid=~Enterotype)
# Here is an example with a different faceting variable.
# Not super useful in this dataset, but always good to explore.
plot_bar(ent10, "SeqTech", fill="Enterotype", facet_grid=~ClinicalStatus)
}
\seealso{
  \code{\link{psmelt}}

  \code{\link{ggplot}}

  \code{\link{qplot}}
}

