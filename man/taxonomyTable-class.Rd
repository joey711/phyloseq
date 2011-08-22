\name{taxonomyTable-class}
\alias{taxonomyTable-class}
\title{An S4 class that holds taxonomic classification data as a character
matrix.}
\description{
  Row indices represent taxa, columns represent taxonomic
  classifiers.
}
\details{
  \describe{ \item{nspecies}{A single positive integer
  specifying the number of distinct species/taxa.
  \code{nspecies} should not be assigned directly by user,
  but is assessed internally during instantiation of
  \code{taxonomyTable} objects.}

  \item{.Data}{This slot is inherited from the
  \code{\link{matrix}} class.} }
}

