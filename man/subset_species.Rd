\name{subset_species}
\alias{subset_species}
\title{Subset species by taxonomic expression}
\usage{
  subset_species(x, ...)
}
\arguments{
  \item{x}{A taxonomyTable class, or \code{phyloseq-class}
  that contains a taxonomyTable. If the \code{taxTab} slot
  is missing in \code{x}, then \code{x} will be returned
  as-is and a warning will be printed to screen.}

  \item{...}{The subsetting expression that should be
  applied to the \code{taxonomyTable}. This is passed on to
  \code{\link{subset}}, and more details and examples about
  how it functions can be found in its documentation.}
}
\value{
  A subsetted object with the same class as \code{x}.
}
\description{
  This is a convenience wrapper around the
  \code{\link{subset}} function. It is intended to speed
  subsetting complex experimental objects with one function
  call. In the case of \code{subset_species}, the
  subsetting will be based on an expression related to the
  columns and values within the \code{taxTab}
  (\code{taxonomyTable} component) slot of \code{x}.
}
\examples{
## ex3 <- subset_species(ex1, Phylum=="Bacteroidetes")
}
\seealso{
  \code{\link{subset_samples}}
}

