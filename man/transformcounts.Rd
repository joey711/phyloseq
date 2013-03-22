\docType{methods}
\name{transform_sample_counts}
\alias{transform_sample_counts}
\alias{transformSampleCounts}
\title{Transform the abundance count data in an \code{otu_table}, sample-by-sample.}
\usage{
  transform_sample_counts(physeq, fun)

  transformSampleCounts(physeq, fun)
}
\arguments{
  \item{physeq}{(Required). \code{\link{phyloseq-class}} of
  \code{\link{otu_table-class}}.}

  \item{fun}{(Required). A single-argument function that
  will be applied to the abundance counts of each sample.
  Can be an anonymous \code{\link[base]{function}}.}
}
\value{
  A transformed \code{otu_table} -- or \code{phyloseq}
  object with its transformed \code{otu_table}. In general,
  trimming is not expected by this method, so it is
  suggested that the user provide only functions that
  return a full-length vector. Filtering/trimming can
  follow, for which the \code{\link{genefilter_sample}} and
  \code{\link{prune_taxa}} functions are suggested.
}
\description{
  This function transforms the sample counts of a taxa
  abundance matrix according to a user-provided function.
  The counts of each sample will be transformed
  individually. No sample-sample interaction/comparison is
  possible by this method.
}
\examples{
#
data(esophagus)
x1 = transform_sample_counts(esophagus, threshrankfun(50))
head(otu_table(x1), 10)
x2 = transform_sample_counts(esophagus, rank)
head(otu_table(x2), 10)
identical(x1, x2)
x3 = otu_table(esophagus) + 5
x3 = transform_sample_counts(x3, log)
head(otu_table(x3), 10)
x4 = transform_sample_counts(esophagus, function(x) round(x^2.2, 0))
head(otu_table(x4), 10)
}
\seealso{
  \code{\link{threshrankfun}}, \code{\link{rank}},
  \code{\link{log}}
}

