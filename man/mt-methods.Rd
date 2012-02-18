\docType{methods}
\name{mt}
\alias{mt}
\alias{mt,otuTable,character-method}
\alias{mt,otuTable,factor-method}
\alias{mt,otuTable,integer-method}
\alias{mt,otuTable,logical-method}
\alias{mt,otuTable,numeric-method}
\alias{mt,phyloseq,ANY-method}
\title{Multiple testing of taxa abundance acccording to sample categories/classes}
\usage{
  mt(physeq, classlabel, minPmaxT="minP", ...)
}
\arguments{
  \item{physeq}{(Required). \code{\link{otuTable-class}} or
  \code{\link{phyloseq-class}}. In this multiple testing
  framework, different taxa correspond to variables
  (hypotheses), and samples to observations.}

  \item{classlabel}{(Required). A single character index of
  the sample-variable in the \code{\link{sampleData}} of
  \code{physeq} that will be used for multiple testing.
  Alternatively, \code{classlabel} can be a custom integer
  (or numeric coercable to an integer), character, or
  factor with length equal to \code{nsamples(physeq)}.

  NOTE: the default test applied to each taxa is a
  two-sample two-sided \code{\link{t.test}}, WHICH WILL
  FAIL with an error if you provide a data variable (or
  custom vector) that contains MORE THAN TWO classes. One
  alternative to consider is an F-test, by specifying
  \code{test="f"} as an additional argument. See the first
  example below, and/or further documentation of
  \code{\link[multtest]{mt.maxT}} or
  \code{\link[multtest]{mt.minP}} for other options and
  formal details.}

  \item{minPmaxT}{(Optional). Character string.
  \code{"mt.minP"} or \code{"mt.maxT"}. Default is to use
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
  Multiple testing of taxa abundance acccording to sample
  categories/classes
}
\examples{
#
## # Simple example, testing genera that sig correlate with Enterotypes
## data(enterotype)
## # Filter samples that don't have Enterotype
## x <- subset_samples(enterotype, !is.na(Enterotype))
## # (the taxa are at the genera level in this dataset)
## mt(x, "Enterotype", test="f")
## # Not surprisingly, Prevotella and Bacteroides top the list.
## # Different test, multiple-adjusted t-test, whether samples are ent-2 or not.
## mt(x, getVariable(x, "Enterotype")==2)
}
\seealso{
  \code{\link[multtest]{mt.maxT}},
  \code{\link[multtest]{mt.minP}}
}

