\docType{methods}
\name{nspecies}
\alias{nspecies}
\alias{nspecies,ANY-method}
\alias{nspecies,otuTable-method}
\alias{nspecies,phylo-method}
\alias{nspecies,phyloseq-method}
\alias{nspecies,taxonomyTable-method}
\title{Get the number of taxa/species.}
\arguments{
  \item{physeq}{\code{\link{phyloseq-class}},
  \code{\link{otuTable-class}},
  \code{\link{taxonomyTable-class}}, or
  \code{\link[ape]{phylo}}}
}
\value{
  An integer indicating the number of taxa / species.
}
\description{
  Get the number of taxa/species.
}
\examples{
#
# # From "picante" package
# data("phylocom")
# tree <- phylocom$phylo
# nspecies(tree)
}
\seealso{
  species.names
}

