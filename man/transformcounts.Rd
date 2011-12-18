\docType{methods}
\name{transformsamplecounts}
\alias{transformsamplecounts}
\alias{transformSampleCounts}
\alias{TransformSampleCounts}
\title{Transform the abundance count data in an \code{otuTable}, sample-by-sample.}
\usage{
  transformsamplecounts(physeq, flist)

  TransformSampleCounts(physeq, flist)

  transformSampleCounts(physeq, flist)
}
\arguments{
  \item{physeq}{(Required). \code{\link{phyloseq-class}} of
  \code{\link{otuTable-class}}.}

  \item{flist}{(Required). A list of functions that will be
  applied to the abundance counts of each sample.}
}
\value{
  A transformed \code{otuTable} -- or \code{phyloseq}
  object with its transformed \code{otuTable}. In general,
  trimming is not expected by this method, so it is
  suggested that the user provide only functions that
  return a full-length vector. Filtering/trimming can
  follow, for which the \code{\link{genefilterSample}} and
  \code{\link{prune_species}} functions are suggested.
}
\description{
  This higher-order function transforms the sample counts
  of a species abundance matrix according to a
  user-provided function or list of functions. The counts
  of each sample will be transformed individually. No
  sample-sample interaction/comparison is possible by this
  method.

  # #

  # #
}
\examples{
#
## data(ex1)
## ex1r <- transformsamplecounts(ex1, rank)
}
\seealso{
  \code{\link{threshrankfun}}, \code{\link{rank}},
  \code{\link{log}}
}

