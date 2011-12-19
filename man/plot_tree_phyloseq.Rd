\name{plot_tree_phyloseq}
\alias{plot_tree_phyloseq}
\title{Plot tree with easy tip annotation.}
\usage{
  plot_tree_phyloseq(physeq, color_factor=NULL,
  shape_factor=NULL, base_size=1, size_scaling_factor =
  0.2, opacity=2/3, custom_color_scale=NULL,
  custom_shape_scale=NULL, type_abundance_value=FALSE,
  printTheseTaxa=NULL, treeTitle="Annotated Tree", ...)
}
\arguments{
  \item{physeq}{(Required). \code{\link{phyloseq-class}}
  with non-empty tree, sampleData, and otuTable
  components.}

  \item{color_factor}{A character string specifying the
  column of the sampleData that will be used for setting
  the color of symbols.}

  \item{shape_factor}{A character string specifying the
  column of the sampleData that will be used for setting
  the shape of symbols.}

  \item{base_size}{The minimum size expansion factor of
  symbols plotted next to tips. The default value is 1.}

  \item{size_scaling_factor}{A numeric, greater than or
  equal to 0, that is multiplied by the log10 of taxa
  abundance; the product of which is summed with the
  \code{base_size} argument to determine the size scaling
  factor provided to \code{\link{tipsymbols}}. The default
  value is 0.15. The larger the value, the larger the
  symbols representing sites with many idividuals of a
  particular taxa. A value of zero means there will be no
  scaling symbol size by the abundance value.}

  \item{opacity}{The opacity (or alpha value). Numeric
  between 0, 1. Defaul value is 2/3.}

  \item{custom_color_scale}{A character vector of the
  desired custom color scale. This should be a scale, not
  an aesthetic map. Therefore, it will in most-cases
  contain only unique elements, unless two different
  categories of data are supposed to have the same color.
  Default value is NULL, which invokes a default color
  scale using the \code{\link{rainbow}} function.}

  \item{custom_shape_scale}{An integer vector of values in
  the categorical scale of symbol shapes, analogous to
  \code{custom_color_scale}. Default value is \code{NULL},
  which uses the fill-able symbols described in
  \code{\link{points}}, beginning with 21.}

  \item{type_abundance_value}{Logical. Whether or not the
  otuTable value (the number of individuals, typically)
  should be added to the center of symbols when the value
  is greater than one. Default is FALSE, indicating no
  labels.}

  \item{printTheseTaxa}{a character vector of the taxa
  names in \code{physeq} that should be labeled on the tree
  plot adjacent to the right. Default is NULL. Not yet
  implemented.}

  \item{treeTitle}{(Optional). Character string, for the
  title of the graphic. Default is \code{"Annotated
  Tree"}.}

  \item{...}{Additional parameters passed on to
  \code{\link{tipsymbols}}.}
}
\value{
  Creates a phylogenetic tree, with additional symbols
  annotated on each tip to indicate in which samples the
  particular taxa was observed.
}
\description{
  Requires a \code{\link{phyloseq-class}} that contains a
  tree (\code{\link{tre}}), sample data
  (\code{\link{sampleData}}), and abundance table
  (\code{\link{otuTable}}).
}
\examples{
# data(ex1)
# ex2 <- ex1
## # for reproducibility
## set.seed(711)
# species.names(ex2) <- sample(species.names(ex1), 50)
# plot_tree_phyloseq(ex2)
# plot_tree_phyloseq(ex2, shape_factor="Diet")
# plot_tree_phyloseq(ex2, color_factor="Gender", shape_factor="Diet")
# plot_tree_phyloseq(ex2, color_factor="Gender", shape_factor="Diet",
# size_scaling_factor=0.6, type_abundance_value=TRUE)
# plot_tree_phyloseq(ex2, color_factor="Gender", shape_factor="Diet",
# size_scaling_factor=0.6, custom_color_scale=c("blue", "magenta"))
# plot(phyloseqTree(ex2), color_factor="Gender", shape_factor="Diet",
# size_scaling_factor=0.6, custom_color_scale=c("blue", "magenta") )
}

