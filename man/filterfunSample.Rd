\name{filterfunSample}
\alias{filterfunSample}
\title{A sample-wise filter function builder, analogous to \code{\link[genefilter]{filterfun}}.}
\usage{
  filterfunSample(...)
}
\arguments{
  \item{...}{A comma-separated list of functions.}
}
\value{
  An enclosure (function) that itself will return a logical
  vector, according to the functions provided in the
  argument list, evaluated in order. The output of
  filterfunSample is appropriate for the `flist' argument
  to the genefilterSample method.
}
\description{
  See the \code{\link[genefilter]{filterfun}}, from the
  Bioconductor repository, for a taxa-/gene-wise filter
  (and further examples).
}
\examples{
## Use simulated abundance matrix
# set.seed(711)
# testOTU <- otuTable(matrix(sample(1:50, 25, replace=TRUE), 5, 5), speciesAreRows=FALSE)
# f1  <- filterfunSample(topk(2))
# wh1 <- genefilterSample(testOTU, f1, A=2)
# wh2 <- c(T, T, T, F, F)
# prune_species(wh1, testOTU)
# prune_species(wh2, testOTU)
}
\seealso{
  \code{\link[genefilter]{filterfun}},
  \code{\link{genefilterSample}}
}

