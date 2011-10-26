\docType{methods}
\name{speciesarerows<-}
\alias{assign-speciesarerows}
\alias{speciesarerows<-}
\alias{speciesarerows<-,otuTable,logical-method}
\alias{speciesarerows<-,phyloseqFather,logical-method}
\title{Manually change speciesAreRows through assignment.}
\arguments{
  \item{x}{An otuTable-class or higher-order object that
  contains an otuTable.}

  \item{value}{A logical of length equal to 1. If
  \code{length(value) > 1}, the additional elements will be
  ignored. Only the first element is assigned to the
  speciesAreRows slot.}
}
\description{
  The speciesAreRows slot is a logical indicating the
  orientation of the abundance table contained in object
  \code{x}.
}
\examples{
#
}

