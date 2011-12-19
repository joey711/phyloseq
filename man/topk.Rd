\name{topk}
\alias{topk}
\title{Make filter fun. the most abundant \code{k} taxa}
\usage{
  topk(k, na.rm=TRUE)
}
\arguments{
  \item{k}{An integer, indicating how many of the most
  abundant taxa should be kept.}

  \item{na.rm}{A logical. Should \code{NA}s be removed.
  Default is \code{TRUE}.}
}
\value{
  Returns a function (enclosure) that will return TRUE for
  each element in the most abundant k values.
}
\description{
  Make filter fun. the most abundant \code{k} taxa
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
  \code{\link{topk}}, \code{\link{topf}},
  \code{\link{topp}}, \code{\link{rm_outlierf}}
}

