\name{sample_sums}
\alias{sample_sums}
\alias{sampleSums}
\title{Returns the total number of individuals observed from each sample.}
\usage{
  sample_sums(x)

  sampleSums(x)
}
\arguments{
  \item{x}{\code{\link{otu_table-class}}, or
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
  but where the orientation of the otu_table is
  automatically handled.
}
\examples{
data(enterotype)
sample_sums(enterotype)
data(esophagus)
sample_sums(esophagus)
}
\seealso{
  \code{\link{taxa_sums}}, \code{\link{rowSums}},
  \code{\link{colSums}}
}

