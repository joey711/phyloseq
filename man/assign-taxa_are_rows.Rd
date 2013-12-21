\docType{methods}
\name{taxa_are_rows<-}
\alias{assign-taxa_are_rows}
\alias{taxa_are_rows<-}
\alias{taxa_are_rows<-,otu_table,logical-method}
\alias{taxa_are_rows<-,phyloseq,logical-method}
\title{Manually change taxa_are_rows through assignment.}
\usage{
  taxa_are_rows(x) <- value
}
\arguments{
  \item{x}{\code{\link{otu_table-class}} or
  \code{\link{phyloseq-class}}}

  \item{value}{A logical of length equal to 1. If
  \code{length(value) > 1}, the additional elements will be
  ignored. Only the first element is assigned to the
  taxa_are_rows slot.}
}
\description{
  The taxa_are_rows slot is a logical indicating the
  orientation of the abundance table contained in object
  \code{x}.
}
\examples{
data(esophagus)
taxa_are_rows(esophagus)
taxa_are_rows(otu_table(esophagus))
}

