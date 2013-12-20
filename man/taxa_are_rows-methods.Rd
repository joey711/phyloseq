\docType{methods}
\name{taxa_are_rows}
\alias{taxa_are_rows}
\alias{taxa_are_rows,ANY-method}
\alias{taxa_are_rows,otu_table-method}
\alias{taxa_are_rows,phyloseq-method}
\title{Access taxa_are_rows slot from otu_table objects.}
\usage{
  taxa_are_rows(physeq)
}
\arguments{
  \item{physeq}{(Required). \code{\link{phyloseq-class}},
  or \code{\link{otu_table-class}}.}
}
\value{
  A logical indicating the orientation of the otu_table.
}
\description{
  Access taxa_are_rows slot from otu_table objects.
}
\seealso{
  \code{\link{otu_table}}
}

