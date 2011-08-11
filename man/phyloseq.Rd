\name{phyloseq}
\alias{phyloseq}
\title{Define the phyloseq class.}
\arguments{
  \item{object}{An object among the set of classes defined
  by the phyloseq package, as well as \code{phylo} and
  \code{phylo4} objects (which are defined by the ape and
  phylobase packages, respectively). For constructing H.O.
  objects, this should be an otuTable object.}

  \item{x1}{(optional) an additional component data object
  (but not otuTable)}

  \item{x2}{(optional) an additional component data object
  (but not otuTable)}

  \item{x3}{(optional) an additional component data object
  (but not otuTable)}
}
\value{
  The class of the returned object depends on the argument
  class(es). To construct a H.O. object, the otuTable must
  be the first object and/or assigned to the named argument
  \code{object}. Otherwise, the order of arguments does not
  matter. If a single component-class object is provided,
  it is simply returned as-is. If a higher-order object is
  provided as input, then the accessor method is invoked,
  and an object of class \code{phyloseq} is returned.
}
\description{
  Contains otuTable and sampleMap slots.

  \code{phyloseq()} is both a constructor and accessor
  method, on account of the ambiguation of `phyloseq' being
  both the package name and a class name (the
  \code{phyloseq} class, which contains an \code{otuTable}
  and \code{sampleMap} slots, only).
}
\details{
  This is the main method suggested for constructing
  \code{phyloseq} higher-order objects from their
  components.
}
\examples{
#
## OTU1 <- otuTable(matrix(sample(0:5,250,TRUE),25,10), speciesAreRows=TRUE)
## tax1 <- taxTab(matrix("abc", 30, 8))
## map1 <- data.frame( matrix(sample(0:3,250,TRUE),25,10),
##   matrix(sample(c("a","b","c"),150,TRUE), 25, 6) )
## map1 <- sampleMap(map1)
## ex1 <- phyloseq(OTU1, map1, tax1)
## phyloseq(ex1)
}
\seealso{
  merge_phyloseq
}
\keyword{constructor}
\keyword{phyloseq}

