\name{sampleData<-}
\alias{assign-sampleData}
\alias{samData<-}
\alias{sampleData<-}
\title{Assign (new) sampleData to \code{x}}
\usage{
  sampleData(x) <- value

  samData(x) <- value
}
\arguments{
  \item{x}{(Required). \code{\link{phyloseq-class}}. The
  object to modify.}

  \item{value}{(Required). Either a
  \code{\link{sampleData-class}}, a \code{data.frame} that
  can be coerced into \code{\link{sampleData-class}}, or a
  \code{\link{phyloseq-class}} that contains a suitable
  \code{sampleData} component to assign to \code{x}. If
  unsure, try \code{\link{sampleData}}\code{(value)}, which
  should return a \code{\link{sampleData-class}} object
  without error.}
}
\value{
  No return. This is an assignment statement.
}
\description{
  This replaces the current \code{sampleData} component of
  \code{x} with \code{value}, if \code{value} is a
  \code{\link{sampleData-class}}. However, if \code{value}
  is a \code{data.frame}, then \code{value} is first
  coerced to a \code{\link{sampleData-class}}, and then
  assigned. Alternatively, if \code{value} is
  \code{\link{phyloseq-class}}, then the
  \code{\link{sampleData}} component will first be accessed
  from \code{value} and then assigned. This makes possible
  some concise assignment/replacement statements when
  adjusting, modifying, or building subsets of
  experiment-level data. See some examples below.
}
\details{
  Internally, this re-builds the
  \code{\link{phyloseq-class}} object using the standard
  \code{\link{phyloseq}} constructor. Thus, index
  mismatches between sample-describing components will not
  be allowed, and subsetting will occurr automatically such
  that only the intersection of sample IDs are included in
  any components. This has the added benefit of re-checking
  (internally) for any other issues.
}
\examples{
#
# data(GlobalPatterns)
# # An example of pruning to just the first 10 samples in GlobalPatterns
# ex2a <- prune_samples(sample.names(GlobalPatterns)[1:10], GlobalPatterns)
# # The following 3 lines produces an ex2b that is equal to ex2a
# ex2b <- GlobalPatterns
# SD <- sampleData(GlobalPatterns)[1:10, ]
# sampleData(ex2b) <- SD
# identical(ex2a, ex2b)
# print(ex2b)
# # Example restoring the original sampleData component. ex2c lacks sampleData
# ex2c <- phyloseq(otuTable(GlobalPatterns), taxTab(GlobalPatterns), tre(GlobalPatterns))
# sampleData(ex2c) <- GlobalPatterns
# identical(ex2c, GlobalPatterns)
# # Can try on ex2b, but other components have only 10 samples. No change.
# sampleData(ex2b) <- GlobalPatterns
# identical(ex2a, ex2b) # still true.
}

