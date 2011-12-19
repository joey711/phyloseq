\name{speciesSums}
\alias{speciesSums}
\title{Returns the total number of individuals observed from each species/taxa/OTU.}
\usage{
  speciesSums(x)
}
\arguments{
  \item{x}{\code{\link{otuTable-class}}, or
  \code{\link{phyloseq-class}}.}
}
\value{
  A \code{\link{numeric-class}} with length equal to the
  number of species in the table, name indicated the taxa
  ID, and value equal to the sum of all individuals
  observed for each taxa in \code{x}.
}
\description{
  A convenience function equivalent to rowSums or colSums,
  but where the orientation of the otuTable is
  automatically handled.
}
\examples{
data(enterotype)
speciesSums(enterotype)
data(esophagus)
speciesSums(esophagus)
}
\seealso{
  \code{\link{sampleSums}}, \code{\link{rowSums}},
  \code{\link{colSums}}
}

