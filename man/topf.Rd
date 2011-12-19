\name{topf}
\alias{topf}
\title{Make filter fun. that returns the top f fraction of taxa in a sample.}
\usage{
  topf(f, na.rm=TRUE)
}
\arguments{
  \item{f}{Single numeric value between 0 and 1.}

  \item{na.rm}{Logical. Should we remove NA values. Default
  \code{TRUE}.}
}
\value{
  A function (enclosure), suitable for
  \code{\link{filterfunSample}}, that will return
  \code{TRUE} for each element in the taxa comprising the
  most abundant f fraction of individuals.
}
\description{
  As opposed to \code{\link{topp}}, which gives the most
  abundant p fraction of observed taxa (richness, instead
  of cumulative abundance. Said another way, topf ensures a
  certain fraction of the total sequences are retained,
  while topp ensures that a certain fraction of
  taxa/species/OTUs are retained.
}
\examples{
# t1 <- 1:10; names(t1)<-paste("t", 1:10, sep="")
# topf(0.6)(t1)
## Use simulated abundance matrix
# set.seed(711)
# testOTU <- otuTable(matrix(sample(1:50, 25, replace=TRUE), 5, 5), speciesAreRows=FALSE)
# f1  <- filterfunSample(topf(0.4))
# (wh1 <- genefilterSample(testOTU, f1, A=1))
# wh2 <- c(T, T, T, F, F)
# prune_species(wh1, testOTU)
# prune_species(wh2, testOTU)
}
\seealso{
  \code{\link{topk}}, \code{\link{topf}},
  \code{\link{topp}}, \code{\link{rm_outlierf}}
}

