\docType{methods}
\name{mt}
\alias{mt}
\alias{mt,otuTable,character-method}
\alias{mt,otuTable,factor-method}
\alias{mt,otuTable,integer-method}
\alias{mt,otuTable,logical-method}
\alias{mt,otuTable,numeric-method}
\alias{mt,phyloseq,ANY-method}
\title{Multiple testing of taxa abundancesa acccording to a sample variate}
\arguments{
  \item{physeq}{(Required). \code{\link{otuTable-class}} or
  \code{\link{phyloseq-class}}.  In this multiple testing
  framework, different taxa correspond to variables
  (hypotheses), and samples to observations.}

  \item{classlabel}{(Required). A single character index of
  the sample-variable in the \code{\link{sampleData}} of
  \code{physeq} that will be used for multiple testing.
  Alternatively, \code{classlabel} can be a custom integer,
  character, or factor with length equal to
  \code{nsamples(physeq)}. In either scenario -- a sample
  variable within the phyloseq-class object or a custom
  vector -- \code{classlabel} must be coercable to a
  logical (i.e. two-level factor). A variable with more
  than two states is not appropriate for this category of
  test.}

  \item{minPmaxT}{(Optional). Character string.
  \code{"mt.minP"} or \code{"mt.maxT"}.  Default is to use
  \code{\link[multtest]{mt.minP}}.}

  \item{...}{(Optional). Additional arguments, forwarded to
  \code{\link[multtest]{mt.maxT}} or
  \code{\link[multtest]{mt.minP}}}
}
\value{
  A dataframe with components specified in the
  documentation for \code{\link[multtest]{mt.maxT}} or
  \code{\link[multtest]{mt.minP}}, respectively.
}
\description{
  Multiple testing of taxa abundancesa acccording to a
  sample variate
}
\examples{
#
## data(ex1)
## mt(ex1, "Diet")
## # Alternatively
## dietfac <- factor(data.frame(sampleData(ex1))[, "Diet"])
## mt(ex1, dietfac)
}
\seealso{
  \code{\link[multtest]{mt.maxT}},
  \code{\link[multtest]{mt.minP}}
}

