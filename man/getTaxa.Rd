\name{getTaxa}
\alias{getTaxa}
\title{Get a unique vector of the observed taxa at a particular taxonomic rank}
\usage{
  getTaxa(physeq, taxonomic.rank=rank.names(physeq)[1])
}
\arguments{
  \item{physeq}{(Required).
  \code{\link{taxonomyTable-class}}, or
  \code{\link{phyloseq-class}}.}

  \item{taxonomic.rank}{(Optional). Character. The
  taxonomic rank to use. Must select from the set indicated
  by \code{getTaxa}. Default is to take the first column of
  the \code{taxonomyTable} component.}
}
\value{
  Character vector. Unique vector of the observed taxa at a
  particular taxonomic rank
}
\description{
  This is a simple accessor function to make it more
  convenient to determine the different taxa present for a
  particular taxonomic rank in a given
  \code{\link{phyloseq-class}} object.
}
\examples{
data(enterotype)
getTaxa(enterotype)
data(GlobalPatterns)
getTaxa(GlobalPatterns, "Family")
}
\seealso{
  getSpecies species.names sample.names getTaxa
}

