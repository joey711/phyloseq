\name{taxTab}
\alias{taxtab}
\alias{taxTab}
\title{Build or access taxonomyTable-class objects.}
\arguments{
  \item{object}{An object among the set of classes defined
  by the phyloseq package that contain taxonomyTable.}
}
\value{
  A taxonomyTable object. It is either grabbed from the
  relevant slot if \code{object} is complex, or built anew
  if \code{object} is a character matrix representing the
  taxonomic classification of species in the experiment.
}
\description{
  \code{taxTab()} is both a constructor and accessor
  method. When the argument is a character matrix, taxTab()
  will attempt to create and return a taxonomyTable-class
  object. In this case, the rows should be named to match
  the \code{species.names} of the other objects to which it
  will ultimately be paired. Alternatively, if the argument
  is an object that contains a taxonomyTable, including a
  taxonomyTable-class object, then the corresponding
  taxonomyTable is returned, as-is.
}
\details{
  This is the main method suggested for constructing
  taxonomyTable objects from a character matrix of taxonomy
  classifications. Each row represents a different species.
  It is also the suggested method for accessing/subsetting
  the taxonomyTable from a more complex object.
  \code{taxTab} is the slot name that holds the
  taxonomyTable-class object in a multi-component phyloseq
  object.
}
\examples{
#
# OTU1 <- otuTable(matrix(sample(0:5,250,TRUE),25,10), speciesAreRows=TRUE)
# tax1 <- taxTab(matrix("abc", 30, 8))
# map1 <- data.frame( matrix(sample(0:3,250,TRUE),25,10),
#   matrix(sample(c("a","b","c"),150,TRUE), 25, 6) )
# map1 <- sampleMap(map1)
# ex1 <- phyloseq(OTU1, map1, tax1)
# taxTab(ex1)
# tax1
}
\seealso{
  otuTable sampleMap tre phyloseq
}
\keyword{OTU}
\keyword{phylo}
\keyword{tre}

