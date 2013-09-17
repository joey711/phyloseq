\name{rarefy_even_depth}
\alias{rarefy_even_depth}
\title{Resample an OTU table such that all samples have the same library size.}
\usage{
  rarefy_even_depth(physeq, sample.size =
  min(sample_sums(physeq)), rngseed = FALSE, replace =
  TRUE, trimOTUs = TRUE)
}
\arguments{
  \item{physeq}{(Required). A \code{\link{phyloseq-class}}
  object that you want to trim/filter.}

  \item{sample.size}{(Optional). A single integer value
  equal to the number of reads being simulated, also known
  as the depth, and also equal to each value returned by
  \code{\link{sample_sums}} on the output.}

  \item{rngseed}{(Optional). A single integer value passed
  to \code{\link{set.seed}}, which is used to fix a seed
  for reproducibly random number generation (in this case,
  reproducibly random subsampling). The default value is
  \code{711}. If set to \code{FALSE}, then no fiddling with
  the RNG seed is performed, and it is up to the user to
  appropriately call \code{\link{set.seed}} beforehand to
  achieve reproducible results.}

  \item{replace}{(Optional). Logical. Whether to sample
  with replacement (\code{TRUE}) or without replacement
  (\code{FALSE}). The default is with replacement
  (\code{replace=TRUE}). Two implications to consider are
  that (1) sampling with replacement is faster and more
  memory efficient as currently implemented; and (2),
  sampling with replacement means that there is a chance
  that the number of reads for a given OTU in a given
  sample could be larger than the original count value, as
  opposed to sampling without replacement where the
  original count value is the maximum possible. Prior to
  phyloseq package version number \code{1.5.20}, this
  parameter did not exist and sampling with replacement was
  the only random subsampling implemented in the
  \code{rarefy_even_depth} function. Note that this default
  behavior was selected for computational efficiency, but
  differs from analogous functions in related packages
  (e.g. subsampling in QIIME).}

  \item{trimOTUs}{(Optional). Logical. Whether to trim OTUs
  from the dataset that are no longer observed in any
  sample (have a count of zero in every sample). The number
  of OTUs trimmed, if any, is printed to standard out as a
  reminder.}
}
\value{
  An object of class \code{phyloseq}. Only the
  \code{otu_table} component is modified.
}
\description{
  This function uses the standard R \code{\link{sample}}
  function to resample from the abundance values in the
  \code{otu_table} component of the first argument,
  \code{physeq}. Often one of the major goals of this
  procedure is to achieve parity in total number of counts
  between samples, as an alternative to other formal
  normalization procedures, which is why a single value for
  the \code{sample.size} is expected. This kind of
  resampling can be performed with and without replacement,
  with replacement being the more
  computationally-efficient, default setting. See the
  \code{replace} parameter documentation for more details.
  We recommended that you explicitly select a random number
  generator seed before invoking this function, or,
  alternatively, that you explicitly provide a single
  positive integer argument as \code{rngseed}.
}
\details{
  This approach is sometimes mistakenly called
  "rarefaction", which
  \href{http://en.wikipedia.org/wiki/Rarefaction}{in
  physics refers to a form of wave decompression;} but in
  this context, ecology, the term refers to a
  \href{http://en.wikipedia.org/wiki/Rarefaction_(ecology)}{repeated
  sampling procedure to assess species richness}, first
  proposed in 1968 by Howard Sanders. In contrast, this
  procedure is not typically repeated (though maybe it
  should be to ensure robust inference), but instead used
  as an \emph{ad hoc} means to normalize microbiome counts
  that have resulted from libraries of widely-differing
  sizes. Here we have intentionally adopted an alternative
  name, \code{rarefy}, that has also been used recently to
  describe this process and, to our knowledge, not
  previously used in ecology.

  Make sure to use \code{\link{set.seed}} for
  exactly-reproducible results of the random subsampling.
}
\examples{
# Test with esophagus dataset
data("esophagus")
esorepT = rarefy_even_depth(esophagus, replace=TRUE)
esorepF = rarefy_even_depth(esophagus, replace=FALSE)
data("GlobalPatterns")
GPrepT = rarefy_even_depth(GlobalPatterns, 1E5, replace=TRUE)
GPrepF = rarefy_even_depth(GlobalPatterns, 1E5, replace=FALSE)
}
\seealso{
  \code{\link{sample}}

  \code{\link{set.seed}}
}

