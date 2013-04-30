\name{psmelt}
\alias{psmelt}
\title{Melt phyloseq data object into large data.frame}
\usage{
  psmelt(physeq)
}
\arguments{
  \item{physeq}{(Required). An
  \code{\link{otu_table-class}} or
  \code{\link{phyloseq-class}}. Function most useful for
  phyloseq-class.}
}
\value{
  A \code{\link{data.frame}}-class table.
}
\description{
  The psmelt function is a specialized melt function for
  melting phyloseq objects (instances of the phyloseq
  class), usually for the purpose of graphics production in
  ggplot2-based phyloseq-generated graphics. It relies
  heavily on the \code{\link[reshape2]{melt}} and
  \code{\link{merge}} functions. Note that ``melted''
  phyloseq data is stored much less efficiently, and so RAM
  storage issues could arise with a smaller dataset
  (smaller number of samples/OTUs/variables) than one might
  otherwise expect. For average-sized datasets, however,
  this should not be a problem. Because the number of OTU
  entries has a large effect on the RAM requirement,
  methods to reduce the number of separate OTU entries, for
  instance by agglomerating based on phylogenetic distance
  using \code{\link{tipglom}}, can help alleviate RAM usage
  problems. This function is made user-accessible for
  flexibility, but is also used extensively by plot
  functions in phyloseq.
}
\examples{
data("GlobalPatterns")
gp.ch = subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
mdf = psmelt(gp.ch)
nrow(mdf)
ncol(mdf)
colnames(mdf)
head(rownames(mdf))
# Create a ggplot similar to
library("ggplot2")
p = ggplot(mdf, aes(x=SampleType, y=Abundance, fill=Genus))
p = p + geom_bar(color="black", stat="identity", position="stack")
print(p)
}
\seealso{
  \code{\link{plot_bar}}

  \code{\link[reshape2]{melt}}

  \code{\link{merge}}
}

