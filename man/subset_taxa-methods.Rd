\docType{methods}
\name{subset_taxa}
\alias{subset_species}
\alias{subset_taxa}
\title{Subset species by taxonomic expression}
\usage{
  subset_taxa(physeq, ...)

  subset_species(physeq, ...)
}
\arguments{
  \item{physeq}{A \code{\link{taxonomyTable-class}}, or
  \code{\link{phyloseq-class}} that contains a
  taxonomyTable. If the \code{tax_table} slot is missing in
  \code{physeq}, then \code{physeq} will be returned as-is
  and a warning will be printed to screen.}

  \item{...}{The subsetting expression that should be
  applied to the \code{taxonomyTable}. This is passed on to
  \code{\link{subset}}, and more details and examples about
  how it functions can be found in its documentation.}
}
\value{
  A subsetted object with the same class as \code{physeq}.
}
\description{
  This is a convenience wrapper around the
  \code{\link{subset}} function. It is intended to speed
  subsetting complex experimental objects with one function
  call. In the case of \code{subset_taxa}, the subsetting
  will be based on an expression related to the columns and
  values within the \code{tax_table} (\code{taxonomyTable}
  component) slot of \code{physeq}.
}
\examples{
## ex3 <- subset_taxa(GlobalPatterns, Phylum=="Bacteroidetes")
}
\seealso{
  \code{\link{subset_samples}}
}

