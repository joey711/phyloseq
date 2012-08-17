\name{estimate_richness}
\alias{estimate_richness}
\title{Summarize richness estimates}
\usage{
  estimate_richness(physeq, split=TRUE)
}
\arguments{
  \item{physeq}{(Required). \code{\link{phyloseq-class}},
  or alternatively, an \code{\link{otu_table-class}}. The
  data about which you want to estimate the richness.}

  \item{split}{(Optional). Logical. Should a separate set
  of richness estimates be performed for each sample? Or
  alternatively, pool all samples and estimate richness of
  the entire set.}
}
\value{
  A \code{data.frame} of the richness estimates, and their
  standard error.
}
\description{
  Performs a number of standard richness estimates, and
  returns the results as a \code{data.frame}. Can operate
  on the cumulative population of all samples in the
  dataset, or by repeating the richness estimates for each
  sample individually. NOTE: You must use untrimmed
  datasets for meaningful results, as these estimates (and
  even the ``observed'' richness) are highly dependent on
  the number of singletons. You can always trim the data
  later on if needed, just not before using this function.
}
\examples{
data(GlobalPatterns)
( S.GP <- estimate_richness(GlobalPatterns) )
# # Make the plots
# plot_richness_estimates(GlobalPatterns, "SampleType")
# plot_richness_estimates(GlobalPatterns, "SampleType", "SampleType")
# For more plotting examples, see plot_richness_estimates()
}
\seealso{
  Check out the custom plotting function,
  \code{\link{plot_richness_estimates}}, for easily showing
  the results of different estimates, with method-specific
  error-bars. Also check out the internal functions
  borrowed from the \code{vegan} package:
  \code{\link[vegan]{estimateR}},
  \code{\link[vegan]{diversity}}
}

