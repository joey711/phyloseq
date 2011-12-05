\name{speciesSums}
\alias{speciesSums}
\title{Returns the total number of individuals observed from each species.}
\usage{
  speciesSums(x)
}
\arguments{
  \item{x}{Any phyloseq object that is or contains an
  otuTable.}
}
\value{
  A named integer vector with length equal to the number of
  species in the table, name indicated the taxa ID, and
  value equal to the sum of all individuals observed for
  each taxa in \code{x}.
}
\description{
  A convenience function equivalent to rowSums or colSums,
  but where the orientation of the otuTable is
  automatically handled. Can take more complex phyloseq
  objects, not just otuTable. Result always derived from
  the abundance values in the associated otuTable, not
  other phyloseq tables
}
\examples{
#
}
\seealso{
  sampleSums rowSums colSums
}

