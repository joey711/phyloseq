\name{taxonomyTable}
\alias{taxonomyTable}
\title{An S4 class that holds taxonomic classification data as a character
matrix.}
\description{
  slots are:
}
\details{
  nspecies -- A single positive integer specifying the
  number of distinct species/taxa. \code{nspecies} should
  not be assigned directly by user, but is assessed
  internally during instantiation of \code{taxonomyTable}
  objects.

  .Data -- This slot is inherited from the
  \code{\link{matrix}} class.
}

