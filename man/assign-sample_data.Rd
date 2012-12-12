\name{sample_data<-}
\alias{assign-sample_data}
\alias{sam_data<-}
\alias{sample_data<-}
\alias{sampleData<-}
\title{Assign (new) sample_data to \code{x}}
\usage{
  sample_data(x) <- value

  sam_data(x) <- value
}
\arguments{
  \item{x}{(Required). \code{\link{phyloseq-class}}. The
  object to modify.}

  \item{value}{(Required). Either a
  \code{\link{sample_data-class}}, a \code{data.frame} that
  can be coerced into \code{\link{sample_data-class}}, or a
  \code{\link{phyloseq-class}} that contains a suitable
  \code{sample_data} component to assign to \code{x}. If
  unsure, try \code{\link{sample_data}}\code{(value)},
  which should return a \code{\link{sample_data-class}}
  object without error.}
}
\value{
  No return. This is an assignment statement.
}
\description{
  This replaces the current \code{sample_data} component of
  \code{x} with \code{value}, if \code{value} is a
  \code{\link{sample_data-class}}. However, if \code{value}
  is a \code{data.frame}, then \code{value} is first
  coerced to a \code{\link{sample_data-class}}, and then
  assigned. Alternatively, if \code{value} is
  \code{\link{phyloseq-class}}, then the
  \code{\link{sample_data}} component will first be
  accessed from \code{value} and then assigned. This makes
  possible some concise assignment/replacement statements
  when adjusting, modifying, or building subsets of
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
# ex2a <- prune_samples(sample_names(GlobalPatterns)[1:10], GlobalPatterns)
# # The following 3 lines produces an ex2b that is equal to ex2a
# ex2b <- GlobalPatterns
# SD <- sample_data(GlobalPatterns)[1:10, ]
# sample_data(ex2b) <- SD
# identical(ex2a, ex2b)
# print(ex2b)
# # Example restoring the original sample_data component. ex2c lacks sample_data
# ex2c <- phyloseq(otu_table(GlobalPatterns), tax_table(GlobalPatterns), phy_tree(GlobalPatterns))
# sample_data(ex2c) <- GlobalPatterns
# identical(ex2c, GlobalPatterns)
# # Can try on ex2b, but other components have only 10 samples. No change.
# sample_data(ex2b) <- GlobalPatterns
# identical(ex2a, ex2b) # still true.
}

