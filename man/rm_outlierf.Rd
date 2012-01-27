\name{rm_outlierf}
\alias{rm_outlierf}
\title{Set to FALSE any outlier species greater than f fractional abundance.}
\usage{
  rm_outlierf(f, na.rm=TRUE)
}
\arguments{
  \item{f}{Single numeric value between 0 and 1. The
  maximum fractional abundance value that a taxa will be
  allowed to have in a sample without being marked for
  trimming.}

  \item{na.rm}{Logical. Should we remove NA values. Default
  \code{TRUE}.}
}
\value{
  A function (enclosure), suitable for
  \code{\link{filterfunSample}}.
}
\description{
  This is for removing overly-abundant outlier taxa, not
  for trimming low-abundance taxa.
}
\examples{
t1 <- 1:10; names(t1)<-paste("t", 1:10, sep="")
rm_outlierf(0.15)(t1)
## Use simulated abundance matrix
# set.seed(711)
# testOTU <- otuTable(matrix(sample(1:50, 25, replace=TRUE), 5, 5), speciesAreRows=FALSE)
# speciesSums(testOTU)
# f1  <- filterfunSample(rm_outlierf(0.1))
# (wh1 <- genefilterSample(testOTU, f1, A=1))
# wh2 <- c(T, T, T, F, F)
# prune_species(wh1, testOTU)
# prune_species(wh2, testOTU)
}
\seealso{
  \code{\link{topk}}, \code{\link{topf}},
  \code{\link{topp}}, \code{\link{rm_outlierf}}
}

