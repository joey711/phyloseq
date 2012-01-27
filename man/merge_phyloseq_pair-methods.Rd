\docType{methods}
\name{merge_phyloseq_pair}
\alias{merge_phyloseq_pair}
\alias{merge_phyloseq_pair,otuTable,otuTable-method}
\alias{merge_phyloseq_pair,phylo,phylo-method}
\alias{merge_phyloseq_pair,sampleData,sampleData-method}
\alias{merge_phyloseq_pair,taxonomyTable,taxonomyTable-method}
\title{Merge pair of phyloseq component data objects of the same class.}
\arguments{
  \item{x}{A character vector of the species in object x
  that you want to keep -- OR alternatively -- a logical
  vector where the kept species are TRUE, and length is
  equal to the number of species in object x. If
  \code{species} is a named logical, the species retained
  is based on those names. Make sure they are compatible
  with the \code{species.names} of the object you are
  modifying (\code{x}).}

  \item{y}{Any \code{phyloseq} object.}
}
\value{
  A single component data object that matches \code{x} and
  \code{y} arguments. The returned object will contain the
  union of the species and/or samples of each. If there is
  redundant information between a pair of arguments of the
  same class, the values in \code{x} are used by default.
  Abundance values are summed for \code{otuTable} objects
  for those elements that describe the same species and
  sample in \code{x} and \code{y}.
}
\description{
  Internal S4 methods to combine pairs of objects of
  classes specified in the phyloseq package. These objects
  must be component data of the same type (class). This is
  mainly an internal method, provided to illustrate how
  merging is performed by the more general
  \code{\link{merge_phyloseq}} function.
}
\details{
  The \code{\link{merge_phyloseq}} function is recommended
  in general.

  Special note: trees are merged using
  \code{\link[ape]{consensus}}.
}
\examples{
#
## # merge two simulated otuTable objects.
## x  <- otuTable(matrix(sample(0:5,200,TRUE),20,10), speciesAreRows=TRUE)
## y  <- otuTable(matrix(sample(0:5,300,TRUE),30,10), speciesAreRows=FALSE)
## xy <- merge_phyloseq_pair(x, y)
## yx <- merge_phyloseq_pair(y, x)
## # merge two simulated taxTab objects
## x <- taxTab(matrix("abc", 20, 6))
## y <- taxTab(matrix("def", 30, 8))
## xy <- merge_phyloseq_pair(x, y)
## # merge two simulated sampleData objects
## x <- data.frame( matrix(sample(0:3,250,TRUE),25,10),
##   matrix(sample(c("a","b","c"),150,TRUE),25,6) )
## x <- sampleData(x)
## y <- data.frame( matrix(sample(4:6,200,TRUE),20,10),
##   matrix(sample(c("d","e","f"),120,TRUE),20,8) )
## y <- sampleData(y)
## merge_phyloseq_pair(x, y)
## data.frame(merge_phyloseq_pair(x, y))
## data.frame(merge_phyloseq_pair(y, x))
}
\seealso{
  \code{\link{merge_phyloseq}} \code{\link{merge_species}}
}

