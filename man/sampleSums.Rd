\name{sampleSums}
\alias{sampleSums}
\title{Returns the total number of individuals observed from each sample.}
\usage{
  sampleSums(x)
}
\arguments{
  \item{x}{\code{\link{otuTable-class}}, or
  \code{\link{phyloseq-class}}.}
}
\value{
  A named \code{\link{numeric-class}} length equal to the
  number of samples in the \code{x}, name indicating the
  sample ID, and value equal to the sum of all individuals
  observed for each sample in \code{x}.
}
\description{
  A convenience function equivalent to rowSums or colSums,
  but where the orientation of the otuTable is
  automatically handled.
}
\examples{
data(enterotype)
sampleSums(enterotype)
data(esophagus)
sampleSums(esophagus)
}
\seealso{
  \code{\link{speciesSums}}, \code{\link{rowSums}},
  \code{\link{colSums}}
}

