\docType{methods}
\name{speciesarerows<-}
\alias{assign-speciesarerows}
\alias{speciesarerows<-}
\alias{speciesarerows<-,otuTable,logical-method}
\alias{speciesarerows<-,phyloseq,logical-method}
\title{Manually change speciesAreRows through assignment.}
\usage{
  speciesarerows(x) <- value
}
\arguments{
  \item{x}{\code{\link{otuTable-class}} or
  \code{\link{phyloseq-class}}}

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
# data(ex1)
# speciesarerows(ex1)
# speciesarerows(otuTable(ex1))
}

