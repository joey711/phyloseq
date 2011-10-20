\name{calcplot}
\alias{calcplot}
\title{Convenience wrapper for performing ordination and plotting.}
\usage{
  calcplot(X, RDA_or_CCA="cca", object=get(all.vars(X)[1]),
  ...)
}
\arguments{
  \item{X}{A formula object, with the left-hand side
  specifying a single phyloseq object that contains (at
  minimum) an \code{otuTable} and a \code{sampleMap}. The
  right-hand side expects at least two different variates,
  present in the \code{sampleMap} of the LHS. For available
  variate names in your object, try
  \code{colnames(sampleMap(ex1))}, where \code{ex1} is the
  phyloseq object containing your data. Because this is a
  formula object, quotes should not be used. See
  \code{\link{formula}} for details about writing a formula
  in R.}

  \item{RDA_or_CCA}{A character string, indicating whether
  the ordination method should be RDA or CCA. Default is
  "cca". Case ignored. Only first letter of string is
  considered.}

  \item{object}{Default is the object specified in the
  left-hand side of \code{X}.}

  \item{...}{Additional plotting arguments, passed on to
  \code{\link{plot_ordination_phyloseq}}}
}
\value{
  A ggplot2 graphics object. If not stored as a variable, a
  graphic object will be produced on the default device.
}
\description{
  Convenience wrapper for performing ordination and
  plotting.
}
\seealso{
  \code{\link{rda.phyloseq}}, \code{\link{cca.phyloseq}},
  \code{\link{rda}}, \code{\link{cca}}
}

