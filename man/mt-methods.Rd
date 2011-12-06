\docType{methods}
\name{mt}
\alias{mt}
\alias{mt,otuTable,factor-method}
\alias{mt,otuTable,integer-method}
\alias{mt,otuTable,numeric-method}
\alias{mt,phyloseq,ANY-method}
\title{Wrapper of mt.maxT and mt.minP for the phyloseq package.}
\usage{
  mt(X, classlabel, minPmaxT="minP", ...)
}
\arguments{
  \item{X}{(Required). An \code{otuTable} or an object that
  contains an \code{otuTable}. In this multiple testing
  framework, different taxa correspond to variables
  (hypotheses), and samples to observations.}

  \item{classlabel}{(Required). Either a single character
  string indicating the column in the \code{sampleMap} of
  X, or, alternatively, an integer vector of length equal
  to the number of samples of X. The latter option is
  equivalent to \code{classlabel} of
  \code{\link[multtest]{mt.maxT}}. A third option is for
  classlabel to be a 2-level factor, with length equal to
  the number of samples of X. If \code{classlabel} is a
  non-integer numeric class, \code{mt()} will attempt to
  coerce it to integer.}

  \item{minPmaxT}{(Optional). Whether to use \code{mt.minP}
  or \code{mt.maxT}. Default is to use \code{mt.minP}}

  \item{...}{(Optional). Additional arguments, forwarded to
  \code{mt.maxT}}
}
\value{
  A dataframe with components specified in the
  documentation for \code{mt.maxT} or \code{mt.minP}
  respectively.
}
\description{
  Wrapper of mt.maxT and mt.minP for the phyloseq package.
}
\examples{
#
## data(ex1)
## mt(ex1, "Diet")
## # Alternatively
## dietfac <- factor(data.frame(sampleMap(ex1))[, "Diet"])
## mt(ex1, dietfac)
}
\seealso{
  \code{\link[multtest]{mt.maxT}},
  \code{\link[multtest]{mt.minP}}
}

