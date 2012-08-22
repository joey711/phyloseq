\docType{methods}
\name{transform_sample_counts}
\alias{transformSampleCounts}
\alias{transform_sample_counts}
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
data(GlobalPatterns)
GP <- GlobalPatterns
## transform_sample_counts can work on phyloseq-class, modifying otu_table only
(GPr <- transform_sample_counts(GP, rank) )
## These two approaches result in identical otu_table
(x1 <- transform_sample_counts( otu_table(GP), threshrankfun(500)) )
(x2 <- otu_table(apply(otu_table(GP), 2, threshrankfun(500)), taxa_are_rows(GP)) )
identical(x1, x2)
}
\seealso{
  \code{\link{threshrankfun}}, \code{\link{rank}},
  \code{\link{log}}
}

