\docType{methods}
\name{rda.phyloseq}
\alias{rda.phyloseq}
\alias{rda.phyloseq,formula-method}
\alias{rda.phyloseq,otuTable-method}
\alias{rda.phyloseq,phyloseqFather-method}
\title{phyloseq package extension for \code{vegan::rda}.}
\usage{
  rda.phyloseq(X, ...)
}
\arguments{
  \item{X}{formula-class object specifying the input. No
  need to splat higher-order objects. \code{rda.phyloseq}
  understands where to find the abundance table and sample
  characteristics. The left-hand side should specify a
  single phyloseq object that contains (at minimum) an
  \code{otuTable} and a \code{sampleMap}. The right-hand
  side expects at least two different variates, present in
  the \code{sampleMap} of the LHS. For available variate
  names in your object, try
  \code{colnames(sampleMap(ex1))}, where \code{ex1} is the
  phyloseq object containing your data. Because this is a
  formula object, quotes should not be used. See
  \code{\link{formula}} for details about writing a formula
  in R. Alternatively, X can be an otuTable-class, or any
  more complex phyloseq-package class that contains an
  otuTable.}

  \item{...}{(Optional). In most cases, if additional
  parameter needs to be specified, it will have the form
  \code{data=DF}, where \code{DF} is a \code{data.frame}
  containing information equivalent to a \code{sampleMap}
  object / component. Only necessary if complex object does
  not already contain \code{sampleMap} or you are keeping
  the data separate for some reason.}
}
\value{
  same output as \code{vegan::\link{rda}}
}
\description{
  A formula is main input to \code{vegan::cca}. This
  complicates dispatch based on object signature. A new
  method with a separate name is defined instead.
}
\examples{
#
## data(ex1)
## rda.phyloseq(ex1 ~ Diet + Gender)
}
\seealso{
  \code{\link{cca.phyloseq}}
}

