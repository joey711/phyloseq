\docType{data}
\name{data-GlobalPatterns}
\alias{data-GlobalPatterns}
\alias{GlobalPatterns}
\title{(Data) Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample (2011)}
\description{
  Published in PNAS in early 2011. This work compared the
  microbial communities from 25 environmental samples and
  three known ``mock communities'' -- a total of 9 sample
  types -- at a depth averaging 3.1 million reads per
  sample. Authors were able to reproduce diversity patterns
  seen in many other published studies, while also
  invesitigating technical issues/bias by applying the same
  techniques to simulated microbial communities of known
  composition.
}
\details{
  abstract from research article (quoted):

  The ongoing revolution in high-throughput sequencing
  continues to democratize the ability of small groups of
  investigators to map the microbial component of the
  biosphere. In particular, the coevolution of new
  sequencing platforms and new software tools allows data
  acquisition and analysis on an unprecedented scale. Here
  we report the next stage in this coevolutionary arms
  race, using the Illumina GAIIx platform to sequence a
  diverse array of 25 environmental samples and three known
  ``mock communities'' at a depth averaging 3.1 million
  reads per sample. We demonstrate excellent consistency in
  taxonomic recovery and recapture diversity patterns that
  were previously reported on the basis of metaanalysis of
  many studies from the literature (notably, the
  saline/nonsaline split in environmental samples and the
  split between host-associated and free-living
  communities). We also demonstrate that 2,000 Illumina
  single-end reads are sufficient to recapture the same
  relationships among samples that we observe with the full
  dataset. The results thus open up the possibility of
  conducting large-scale studies analyzing thousands of
  samples simultaneously to survey microbial communities at
  an unprecedented spatial and temporal resolution.

  (end quote)

  Many thanks to J. Gregory Caporaso for directly providing
  the OTU-clustered data files for inclusion in this
  package.
}
\examples{
# data(GlobalPatterns)
# # Load the GlobalPatterns dataset into the workspace environment
# data(GlobalPatterns)
# # Look at the different values for SampleType
# getVariable(GlobalPatterns, "SampleType")
# ################################################################################
# # Reproduce Figure 4 from the article, but using Jaccard distance,
# # and different clustering methods (UPGMA=="average" used in article)
# # The default method for hclust() uses complete-linkage clustering (method="complete")
# ################################################################################
# # Calculate the jaccard distance between each sample
# jaccdist <- distance(GlobalPatterns, "jaccard")
# plot(hclust(jaccdist, "average"), labels=getVariable(GlobalPatterns, "SampleType"))
# # A different method ("complete-linkage")
# plot(hclust(jaccdist), labels=getVariable(GlobalPatterns, "SampleType"), col=cols)
# # In case you decide to color the tip labels
# colorScale <- rainbow(length(levels(getVariable(GlobalPatterns, "SampleType"))))
# cols       <- colorScale[getVariable(GlobalPatterns, "SampleType")]

# ################################################################################
# # Reproduce Figure 5, but in 2-D
# ################################################################################
# coords <- pcoa(UniFrac(GlobalPatterns))$vectors
# DF     <- data.frame(sampleData(GlobalPatterns), coords)
# ggplot(DF, aes(x=Axis.1, y=Axis.2, color=SampleType)) +
# geom_point(size=4) +
# geom_line() +
# opts(title = "PCoA on unweighted UniFrac distance")

# ################################################################################
# # Reproduce Figure 5 (but in 2-D and using jaccard distance / nmMDS)
# ################################################################################
# # Choose number of axes for non-metric MDS
# N <- 2
# # Perform non-metric multi-dimensional scaling, 3 axes (k=3)
# coords <- scores(metaMDS(jaccdist, k=N))
# # Add the NMDS coordinates to the sample data.frame, DF
# DF    <- data.frame(sampleData(GlobalPatterns), coords)
# # plot the MDS of jaccard-distances, and shade points by soil treatments
# # (two axes only, 3-axes used in Fig 5)
# ggplot(DF, aes(x=NMDS1, y=NMDS2, color=SampleType)) +
# geom_point(size=4) +
# geom_line() +
# opts(title = paste("nmMDS on Jaccard distance, ", N, " axes", sep=""))

# ################################################################################
# # Reproduce Figure 5 (but use Jaccard distance / PCoA)
# ################################################################################
# # use principle coordinates analysis (as in article)
# coords <- pcoa(jaccdist)$vectors

# # Add the PCoA coordinates to the sample data.frame, DF
# DF    <- data.frame(sampleData(GlobalPatterns), coords)

# # plot the PCoA on jaccard-distances, and shade points by soil treatments
# # (First-two axes only, could show 3 as in Fig 5, if desired)
# ggplot(DF, aes(x=Axis.1, y=Axis.2, color=SampleType)) +
# geom_point(size=4) +
# geom_line() +
# opts(title = paste("PCoA on Jaccard distance, two axes", sep=""))

# ################################################################################
# # Reproduce Figure 5, but using correspondence analysis
# ################################################################################
# gpdca  <- ordinate(GlobalPatterns, "DCA")
# coords <- scores(gpdca)$sites
# DF     <- data.frame(sampleData(GlobalPatterns), coords)
# ggplot(DF, aes(x=CA1, y=CA2, color=SampleType)) +
# geom_point(size=4) +
# geom_line() +
# opts(title = paste("DCA on abundances, first two axes", sep=""))
}
\author{
  Caporaso, J. G., et al.
}
\references{
  Caporaso, J. G., et al. (2011). Global patterns of 16S
  rRNA diversity at a depth of millions of sequences per
  sample. PNAS, 108, 4516-4522. PMCID: PMC3063599

  The primary article can be viewed/downloaded at:
  \url{http://www.pnas.org/content/108/suppl.1/4516.short}
}
\keyword{data}

