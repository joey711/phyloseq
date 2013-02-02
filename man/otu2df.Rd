\name{otu2df}
\alias{otu2df}
\title{DEPRECATED. SEE \code{\link{psmelt}} converts OTU-table to data.frame}
\usage{
  otu2df(physeq, taxavec, map, keepOnlyTheseTaxa=NULL,
  threshold=NULL)
}
\arguments{
  \item{physeq}{An \code{otu_table} object.}

  \item{taxavec}{A character vector of the desired
  taxonomic names to categorize each species in physeq.}

  \item{map}{The corresponding sample_data object for
  \code{physeq}.}

  \item{keepOnlyTheseTaxa}{A vector of the taxonomic labels
  that you want included. If NULL, the default, then all
  taxonomic labels are used, except for the empty character
  string, ``'', which is trimmed away.}

  \item{threshold}{A [0,1] numeric. Fraction of abundance
  of the taxonomic groups to keep for each sample. The
  higher the value, the larger the diversity of taxonomic
  groups included. That is, a greater number of the rare
  groups are included. If NULL (or 1), the default, all
  taxonomic groups are included.}
}
\description{
  Was used for plotting with \code{\link{plot_taxa_bar}} in
  the ggplot2 framework.
}
\keyword{internal}

