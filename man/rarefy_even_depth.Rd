\name{rarefy_even_depth}
\alias{rarefy_even_depth}
\title{Perform a random subsampling of an OTU table to a level of even depth.}
\usage{
  rarefy_even_depth(physeq, sample.size =
  min(sample_sums(physeq)), rngseed = 711)
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
}
\value{
  An object of class \code{phyloseq}. Only the
  \code{otu_table} component is modified.
}
\description{
  This function uses the \code{\link{sample}} function to
  randomly subset from the abundance values in each sample
  of the \code{otu_table} component in the \code{physeq}
  argument. Sampling is performed with replacement from a
  vector of taxa indices, with length equal to the argument
  to \code{sample.size}, and probability according to the
  abundances for that sample in \code{physeq}.
}
\details{
  This approach is sometimes (somewhat mistakenly) called
  "rarefaction", or "rarefying", but it is actually a
  single random sample-wise subsampling procedure. The
  original rarefaction procedure includes many random
  subsampling iterations at increasing depth as a means to
  infer richness/alpha-diversity based on apparent
  convergence.

  Make sure to use \code{\link{set.seed}} for
  exactly-reproducible results of the random subsampling.
}
\examples{
# Test with esophagus dataset
data("esophagus")
eso  = rarefy_even_depth(esophagus)
eso0 = prune_taxa(taxa_names(eso), esophagus)
plot(as(otu_table(eso), "vector"), as(otu_table(eso0), "vector"))
UniFrac(eso); UniFrac(eso0); UniFrac(esophagus)
# Test with GlobalPatterns dataset
data("GlobalPatterns")
GP.chl <- subset_taxa(GlobalPatterns, Phylum=="Chlamydiae")
# remove the samples that have less than 20 total reads from Chlamydiae
GP.chl <- prune_samples(names(which(sample_sums(GP.chl)>=20)), GP.chl)
# # (p <- plot_tree(GP.chl, color="SampleType", shape="Family", label.tips="Genus", size="abundance"))
GP.chl.r = rarefy_even_depth(GP.chl)
GP.chl = prune_taxa(taxa_names(GP.chl.r), GP.chl)
plot(as(otu_table(GP.chl.r), "vector"), as(otu_table(GP.chl), "vector"))
# Try ordination of GP.chl and GP.chl.r (default distance is unweighted UniFrac)
plot_ordination(GP.chl, ordinate(GP.chl, "MDS"), color="SampleType") #+ geom_point(size=5)
plot_ordination(GP.chl.r, ordinate(GP.chl.r, "MDS"), color="SampleType") #+ geom_point(size=5)
}
\seealso{
  \code{\link{sample}}

  \code{\link{set.seed}}
}

