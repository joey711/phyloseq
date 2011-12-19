\name{calcplot}
\alias{calcplot}
\title{Convenience wrapper for performing ordination and plotting.}
\usage{
  calcplot(X, RDA_or_CCA="cca", object=get(all.vars(X)[1]),
  ...)
}
\arguments{
  \item{X}{(Required). A formula object. The left-hand side
  specifying a single phyloseq object that contains (at
  minimum) an \code{otuTable} and a \code{sampleData}. The
  right-hand side should contain the label of at least one
  variate present in the \code{sampleData} of the LHS
  object. There is one supported alternative, the special
  \code{"."} for the RHS, in which case all available
  variates will be specified as constraints in the
  constrained ordination, but not shaded. This behavior is
  less explicit and may depend somewhat arbitrarily on the
  order of variates in your \code{sampelMap}; as such, it
  should be used with caution.

  For example, there are two variates in the example object
  obtained with \code{data(ex1)}. \code{X} can be specified
  as

  \code{ex1 ~ Diet + Gender}

  for using both as constraints in the ordination (and
  shading), or as

  \code{ex1 ~ Diet}

  if you only wanted to constrain on the Diet variable, for
  example.

  For available variate names in your object, try
  \code{colnames(sampleData(ex1))}, where \code{ex1} is the
  phyloseq object containing your data. Because this is a
  formula object, quotes should not be used. See
  \code{\link{formula}} for details about writing a formula
  in \code{R}.}

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
\examples{
## data(ex1)
## calcplot(ex1 ~ Diet + Gender)
}
\seealso{
  \code{\link{rda.phyloseq}}, \code{\link{cca.phyloseq}},
  \code{\link{rda}}, \code{\link{cca}}
}

