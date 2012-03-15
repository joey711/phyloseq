\name{plot_taxa_bar}
\alias{plot_taxa_bar}
\alias{taxaplot}
\title{Create a structured barplot graphic of the taxonomic groups.}
\usage{
  plot_taxa_bar(otu, taxavec="Domain",
  showOnlyTheseTaxa=NULL, threshold=NULL,
  x_category="sample", fill_category=x_category,
  facet_formula = . ~ TaxaGroup, OTUpoints=FALSE,
  labelOTUs=FALSE)

  taxaplot(otu, taxavec = "Domain", showOnlyTheseTaxa =
  NULL, threshold = NULL, x_category = "sample",
  fill_category = x_category, facet_formula = . ~
  TaxaGroup, OTUpoints = FALSE, labelOTUs = FALSE)
}
\arguments{
  \item{otu}{(Required). An \code{otuTable} object, or
  higher-order object that contains an otuTable and
  sampleData (e.g. ``otuSam'' class and its superclasses.).
  If \code{otu} does not contain a taxTab slot (is a class
  that does not have ``Tax'' in its title), then the second
  argument, \code{taxavec}, is required and should have
  length equal to the number of species/taxa in
  \code{otu}.}

  \item{taxavec}{A character vector of the desired
  taxonomic names to categorize each species in \code{otu}.
  If \code{otu} is a higher-order object that contains a
  taxonomyTable, then taxavec can alternatively specify the
  desired taxonomic level as a character string of length
  1. E.g. \code{taxavec = "Phylum"}. Default value is
  \code{"Domain"}.}

  \item{showOnlyTheseTaxa}{A vector of the taxonomic labels
  that you want included. If NULL, the default, then all
  taxonomic labels are used, except for the empty character
  string, ``'', which is trimmed away.}

  \item{threshold}{A [0,1] numeric. Fraction of abundance
  of the taxonomic groups to keep for each sample. The
  higher the value, the larger the diversity of taxonomica
  groups included. That is, a greater number of the rare
  groups are included. If NULL (or 1), the default, all
  taxonomic groups are included.}

  \item{x_category}{A character string indicating which
  sampleData column should be used to define the horizontal
  axis categories. Default is \code{"sample"}. Note that a
  few column-names are added by default and are available
  as options. They are ``sample'', ``Abundance'', and
  ``TaxaGroup''.}

  \item{fill_category}{A character string indicating which
  sampleData column should be used to define the fill color
  of the bars. This does not have to match
  \code{x_category}, but does so by default. Note that a
  few column-names are added by default and are available
  as options. They are ``sample'', ``Abundance'', and
  ``TaxaGroup''.}

  \item{facet_formula}{A formula object as used by
  \code{\link{facet_grid}} in \code{\link{ggplot}} or
  \code{\link{qplot}} commands The default is: \code{. ~
  TaxaGroup}. Note that a few column-names are added by
  default and are available as options. They are
  ``sample'', ``Abundance'', and ``TaxaGroup''. E.g. An
  alternative \code{facet_grid} could be \code{sample ~
  TaxaGroup}.}

  \item{OTUpoints}{(Optional). Logical. Default
  \code{FALSE}. Whether to add small grey semi-transparent
  points for each OTU. Helps convey the relative
  distribution within each bar if it combines many
  different OTUs. For datasets with large numbers of
  samples and for complicated plotting arrangements, this
  might be too cluttered to be meaningful.}

  \item{labelOTUs}{(Optional). Logical. Default
  \code{FALSE}. Whether to add a label over the top few
  OTUs within each bar. As with \code{OTUpoints}, this is
  probably not a good idea for plots with large complexity.
  For low numbers of total OTUs this can be informative,
  and help display multiple layers of information on the
  same graphic.}
}
\value{
  A ggplot2 graphic object.
}
\description{
  This function wraps \code{ggplot2} plotting, and returns
  a \code{ggplot2} graphic object that can be saved or
  further modified with additional layers, options, etc.
  The main purpose of this function is to quickly and
  easily create informative summary graphics of the
  differences in taxa abundance between samples in an
  experiment.
}
\details{
  The vertical axis is always relative abundance, but the
  data can be further organized at the horizontal axis and
  faceting grid by any combination of variates present in
  the sampleData component of \code{otu}.
}
\examples{
##
# data(enterotype)
# TopNOTUs <- names(sort(speciesSums(enterotype), TRUE)[1:10])
# ent10   <- prune_species(TopNOTUs, enterotype)
# (p <- plot_taxa_bar(ent10, "Genus", x="SeqTech", fill="TaxaGroup") +
#    facet_wrap(~Enterotype) )
}
\seealso{
  \code{\link{otu2df}}, \code{\link{qplot}},
  \code{\link{ggplot}}
}

