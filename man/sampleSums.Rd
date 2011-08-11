\name{sampleSums}
\alias{sampleSums}
\title{Returns the total number of individuals observed from each sample.}
\arguments{
  \item{x}{Any phyloseq-package object that is or contains
  an otuTable.}
}
\value{
  the total number of individuals present in each sample. A
  named integer vector with length equal to the number of
  samples in the table, name indicated the sample ID, and
  value equal to the sum of all individuals observed for
  each sample in \code{x}.
}
\description{
  A convenience function equivalent to rowSums or colSums,
  but where the orientation of the otuTable is
  automatically handled. Can take more complex phyloseq
  objects, not just otuTable. Result always derived from
  the abundance values in the associated otuTable, not
  other phyloseq tables.
}
\examples{
#
}
\seealso{
  speciesSums rowSums colSums sum
}
\keyword{samples}
\keyword{sum}
\keyword{sums}

