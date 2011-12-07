\name{otuTable-class}
\alias{otuTable-class}
\title{The S4 class for storing taxa-abundance information.}
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
  \describe{ \item{speciesAreRows}{ A single logical
  specifying the orientation of the abundance table. }

  \item{.Data}{This slot is inherited from the
  \code{\link{matrix}} class.} }
}

