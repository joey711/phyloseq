\docType{methods}
\name{transformsamplecounts}
\alias{transformsamplecounts}
\alias{transformSampleCounts}
\alias{TransformSampleCounts}
\title{Transform the abundance count data in an \code{otuTable}, sample-by-sample.}
\usage{
  transformsamplecounts(physeq, fun)

  TransformSampleCounts(physeq, fun)

  transformSampleCounts(physeq, fun)
}
\arguments{
  \item{physeq}{(Required). \code{\link{phyloseq-class}} of
  \code{\link{otuTable-class}}.}

  \item{fun}{(Required). A single-argument function that
  will be applied to the abundance counts of each sample.
  Can be an anonymous \code{\link[base]{function}}.}
}
\value{
  A transformed \code{otuTable} -- or \code{phyloseq}
  object with its transformed \code{otuTable}.  In general,
  trimming is not expected by this method, so it is
  suggested that the user provide only functions that
  return a full-length vector. Filtering/trimming can
  follow, for which the \code{\link{genefilterSample}} and
  \code{\link{prune_species}} functions are suggested.
}
\description{
  This function transforms the sample counts of a species
  abundance matrix according to a user-provided function.
  The counts of each sample will be transformed
  individually. No sample-sample interaction/comparison is
  possible by this method.

  # #

  # #
}
\examples{
#
data(ex1)
## transformsamplecounts can work on phyloseq-class, modifying otuTable only
(ex1r <- transformsamplecounts(ex1, rank) )
## These two approaches result in identical otuTable
(x1 <- transformsamplecounts( otuTable(ex1), threshrankfun(500)) )
(x2 <- otuTable(apply(otuTable(ex1), 2, threshrankfun(500)), speciesAreRows(ex1)) )
identical(x1, x2)
}
\seealso{
  \code{\link{threshrankfun}}, \code{\link{rank}},
  \code{\link{log}}
}

