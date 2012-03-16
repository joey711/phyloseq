\name{filter_taxa}
\alias{filter_taxa}
\title{Filter taxa based on abundance criteria}
\usage{
  filter_taxa(physeq, flist, prune=FALSE)
}
\arguments{
  \item{physeq}{(Required). A \code{\link{phyloseq-class}}
  object that you want to trim/filter.}

  \item{flist}{(Required). A function or list of functions
  that take a vector of abundance values and return a
  logical. Some canned useful function types are included
  in the \code{genefilter}-package.}

  \item{prune}{(Optional). A logical. Default \code{FALSE}.
  If \code{TRUE}, then the function returns the pruned
  \code{\link{phyloseq-class}} object, rather than the
  logical vector of taxa that passed the filter.}
}
\value{
  A logical vector equal to the number of species (taxa) in
  \code{physeq}. This can be provided directly to
  \code{\link{prune_species}} as first argument.
  Alternatively, if \code{prune==TRUE}, the pruned
  \code{\link{phyloseq-class}} object is returned instead.
}
\description{
  This is analogous to \code{\link[genefilter]{genefilter}}
  for microarray filtering. Basically an extension of
  \code{\link[genefilter]{genefilter}} (from the
  Bioconductor repository) for phyloseq objects.
}
\examples{
# library("genefilter")
# data("enterotype")
# flist <- filterfun(kOverA(5, 2e-08), allNA)
# ans   <- filter_taxa(enterotype, flist)
# trimmed.enterotype <- prune_species(ans, enterotype)
# sum(!ans); nspecies(trimmed.enterotype)
# filter_taxa(enterotype, flist, TRUE)
}
\seealso{
  \code{\link[genefilter]{filterfun}},
  \code{\link{genefilterSample}},
  \code{\link{filterfunSample}}
}

