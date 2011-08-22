\name{otuTable-class}
\alias{otuTable-class}
\title{The S4 class that holds taxa-abundance information.}
\description{
  Because orientation of these tables can vary by method,
  the orientation is defined explicitly in the
  \code{speciesAreRows} slot (a logical). The
  \code{otuTable} class inherits the \code{\link{matrix}}
  class to store abundance values. Various standard subset
  and assignment nomenclature has been extended to apply to
  the \code{otuTable} class, including square-bracket,
  \code{\link{t}}, etc.
}
\details{
  \describe{ \item{speciesAreRows}{A single logical
  specifying the orientation of the abundance table.}

  \item{nspecies}{A single positive integer specifying the
  number of distinct species/taxa.}

  \item{nsamples}{ A single positive integer specifying the
  number of distinct samples. Both \code{nsamples} and
  \code{nspecies} should not be assigned directly by user,
  but are assessed internally during initialization of
  \code{otuTable} objects.}

  \item{species.names}{ A character vector of the names of
  the taxa / species, so that they can be accessed,
  modified, subsetted, without the need for the user to
  consider the orientation of the abundance table.}

  \item{sample.names}{A character vector of the names of
  the samples, so that they can be accessed, modified,
  subsetted, without the need for the user to consider the
  orientation of the abundance table. Both
  \code{sample.names} and \code{species.names} can be
  assigned, e.g. \code{species.names(object)<-}, however,
  this has the effect of internally replacing the otuTable
  to be a subset that matches the new name vector, rather
  than replacing the index names in-place.}

  \item{.Data}{This slot is inherited from the
  \code{\link{matrix}} class.} }
}

