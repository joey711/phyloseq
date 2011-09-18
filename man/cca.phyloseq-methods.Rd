\docType{methods}
\name{cca.phyloseq}
\alias{cca.phyloseq}
\alias{cca.phyloseq,formula-method}
\alias{cca.phyloseq,otuTable-method}
\alias{cca.phyloseq,phyloseqFather-method}
\title{phyloseq package extension for \code{vegan::cca}.}
\arguments{
  \item{X}{formula-class object specifying the input. No
  need to splat higher-order objects. \code{cca.phyloseq}
  understands where to find the abundance table and sample
  characteristics. Alternatively, X can be an otuTable, or
  any more complex phyloseq-package class that contains an
  otuTable.}

  \item{data}{\code{data.frame} containing information
  equivalent to a \code{sampleMap} object / component. Only
  necessary if complex object does not already contain
  \code{sampleMap} or you are keeping the data separate for
  some reason.}
}
\value{
  same output as \code{vegan::cca}.
}
\description{
  A formula is main input to \code{vegan::cca}. This
  complicates dispatch based on object signature. A new
  method with a separate name is defined instead.
}
\examples{
#
## data(ex1)
## cca.phyloseq(ex1 ~ Diet + Gender)
}
\seealso{
  \code{\link{rda.phyloseq}}
}

