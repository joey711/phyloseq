\docType{methods}
\name{transformsamplecounts}
\alias{transformsamplecounts}
\alias{transformSampleCounts}
\alias{TransformSampleCounts}
\title{Transform the abundance count data in an \code{otuTable}, sample-by-sample.}
\usage{
  transformsamplecounts(x, flist)

  TransformSampleCounts(x, flist)

  transformSampleCounts(x, flist)
}
\arguments{
  \item{x}{the \code{otuTable} or higher-order object that
  contains an \code{otutable}.}

  \item{flist}{list of functions that will be applied to
  the abundance counts of each sample.}
}
\value{
  The transformed \code{otuTable} or higher-order object
  with its \code{otuTable} transformed. In general,
  trimming is not expected by this method, so it is
  suggested that the user provide only functions that
  return a full-length vector.
}
\description{
  This higher-order function transforms the sample counts
  of a species abundance matrix according to a
  user-provided function or list of functions. The counts
  of each sample will be transformed individually. No
  sample-sample interaction/comparison is possible by this
  method.
}
\examples{
#
## data(ex1)
## ex1r <- transformsamplecounts(ex1, rank)
}

