\name{rank.names}
\alias{rank.names}
\title{Get the names of the taxonomic ranks}
\usage{
  rank.names(physeq, errorIfNULL=TRUE)
}
\arguments{
  \item{physeq}{(Required).
  \code{\link{taxonomyTable-class}}, or
  \code{\link{phyloseq-class}}.}

  \item{errorIfNULL}{(Optional). Logical. Should the
  accessor stop with an error if the slot is empty
  (\code{NULL})? Default \code{TRUE}.}
}
\value{
  Character vector. The names of the available taxonomic
  ranks.
}
\description{
  This is a simple accessor function to make it more
  convenient to determine the taxonomic ranks that are
  available in a given \code{\link{phyloseq-class}} object.
}
\examples{
data(enterotype)
rank.names(enterotype)
}
\seealso{
  getSpecies species.names sample.names getTaxa
}

