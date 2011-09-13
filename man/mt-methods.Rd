\docType{methods}
\name{mt}
\alias{mt}
\alias{mt,otuSam,character-method}
\alias{mt,otuSam,factor-method}
\alias{mt,otuSam,integer-method}
\alias{mt,otuSamTax,ANY-method}
\alias{mt,otuSamTree,ANY-method}
\alias{mt,otuTable,ANY-method}
\alias{mt,otuTable,factor-method}
\alias{mt,otuTree,ANY-method}
\title{Wrapper of mt.maxT and mt.minP for the phyloseq package.}
\arguments{
  \item{X}{An \code{otuTable} or an object that contains an
  \code{otuTable}. In this multiple testing framework,
  different taxa correspond to variables (hypotheses), and
  samples to observations.}

  \item{classlabel}{Either a single character string
  indicating the column in the \code{sampleMap} of X, or,
  alternatively, an integer vector of length equal to the
  number of samples of X. The latter option is equivalent
  to \code{classlabel} of \code{mt.maxT} in
  \code{multtest}. A third option is for classlabel to be a
  2-level factor, with length equal to the number of
  samples of X.}

  \item{...}{Additional arguments, forwarded to
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
## mt(otuSamTax(ex1), "Diet", minPmaxT="minsldzxc")
}
\seealso{
  mt.maxT mt.minP
}

