\name{plot_ordination_phyloseq}
\alias{plot_ordination_phyloseq}
\title{Convenient rendering of ordination results using ggplot2.}
\usage{
  plot_ordination_phyloseq(mod, object, plot_title =
  as(mod$call, "character")[1], man.colors=NULL,
  man.sizes=NULL, man.shapes=NULL, species_alpha=1/10,
  sites_alpha=2/3, species_color_category =NULL,
  species_shape_category =NULL, species_size_category
  =NULL, site_color_category=NULL,
  site_shape_category=NULL, site_size_category=NULL,
  add_sites_data = NULL, add_taxa_data = NULL)
}
\arguments{
  \item{mod}{A \code{cca} or \code{rda} results object. See
  \code{\link{cca.object}}. For phyloseq objects, you
  probably want to see \code{\link{cca.phyloseq}} or
  \code{\link{rda.phyloseq}}.}

  \item{object}{A \code{phyloseq} object or one of its
  child classes (superclasses). Required. Necessary because
  ordination result objects will only refer to the
  abundance table, not the phyloseq object used to create
  them.}

  \item{plot_title}{Default is to use the name of the
  ordination method, stored in the ordination result
  object, accessed via \code{as(mod$call, "character")[1]}.
  Any character string will suffice.}

  \item{man.colors}{A character vector of colors to overide
  the default color scale used by ggplot2. Default is
  \code{NULL}.}

  \item{man.sizes}{An integer vector of sizes to overide
  the default size scale used by ggplot2. Default is
  \code{NULL}.}

  \item{man.shapes}{An integer vector of shape indices to
  overide the defaults used by ggplot2. Default is
  \code{NULL}.}

  \item{species_alpha}{A single numeric. The alpha
  (transparency) value to use for all \code{geom_points}
  that represent taxa. Default is \code{1/10}.}

  \item{sites_alpha}{A single numeric. The alpha
  (transparency) value to use for all \code{geom_points}
  that represent samples/sites. Default is \code{2/3}.}

  \item{species_color_category}{A manual color scale for
  taxa. Default is \code{NULL}.}

  \item{species_shape_category}{A manual shape scale for
  taxa. Default is \code{NULL}.}

  \item{species_size_category}{A manual size scale for
  taxa. Default is \code{NULL}.}

  \item{site_color_category}{A manual color scale for
  samples/sites. Default is \code{NULL}.}

  \item{site_shape_category}{A manual shape scale for
  samples/sites. Default is \code{NULL}.}

  \item{site_size_category}{A manual size scale for
  samples/sites. Default is \code{NULL}.}

  \item{add_sites_data}{A \code{data.frame} object
  providing additional data about the samples/sites that is
  not already available in the \code{sampleMap} of your
  phyloseq object specified by \code{object}. As such, it
  must have the same number of rows as as there are samples
  in \code{object}. Default is \code{NULL}.}

  \item{add_taxa_data}{A \code{matrix} object providing
  additional data about the taxa that is not already
  available in the \code{taxTab} of your phyloseq object
  specified by \code{object}, if it has a \code{taxTab}
  slot. It is not required. As such, it must have the same
  number of rows as there are taxa in \code{object}.
  Default is \code{NULL}.}
}
\value{
  A ggplot2 plot object.
}
\description{
  This convenience function includes many useful defaults
  for generating a plot of your ordination. This
  supplements the default plots created by the plot method
  extension from the vegan package, \code{\link{plot.cca}}.
  For further examples, see the phyloseq vignette.
}
\details{
  Because the ggplot2 objects can be stored as variables,
  and then further modified, it is not necessary to support
  additional arguments as \code{...}. If further
  modifications to the plot produced by this function are
  desired, they can be provided as parameter updates to the
  returned object. See the specific examples provided by
  the phyloseq vignette, or the more comprehensive
  documentation available for the ggplot2 package itself.
}
\seealso{
  plot.cca cca.object
}

