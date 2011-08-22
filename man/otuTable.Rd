\name{otuTable}
\alias{otutable}
\alias{otuTable}
\title{Build or access otuTable objects.}
\arguments{
  \item{object}{An object among the set of classes defined
  by the phyloseq package that contain an otuTable.}

  \item{...}{(optional) ignored unless object is a matrix.
  speciesAreRows is specified here. Must be a named
  argument.}
}
\value{
  An otuTable object. It is either grabbed from the
  relevant slot if \code{object} is complex, or built anew
  if \code{object} is an integer matrix representing the
  species-abundance table.
}
\description{
  \code{otuTable()} is both a constructor and accessor
  method. When the first argument is a matrix, otuTable()
  will attempt to create and return an otuTable-class
  object, which further depends on whether or not
  speciesAreRows is provided as an additional argument.
  Alternatively, if the first argument is an object that
  contains an otuTable, including an otuTable-class object,
  then the corresponding otuTable is returned, as-is.
}
\details{
  This is the main method suggested for constructing
  otuTable objects from an abundance matrix. It is also the
  suggested method for accessing subsetting the otuTable
  from a more complex object.
}
\examples{
#
# OTU1 <- otuTable(matrix(sample(0:5,250,TRUE),25,10), speciesAreRows=TRUE)
# tax1 <- taxTab(matrix("abc", 30, 8))
# map1 <- data.frame( matrix(sample(0:3,250,TRUE),25,10),
#    matrix(sample(c("a","b","c"),150,TRUE), 25, 6) )
# map1 <- sampleMap(map1)
# ex1 <- phyloseq(OTU1, map1, tax1)
# otuTable(ex1)
}
\seealso{
  sampleMap taxTab tre phyloseq
}
\keyword{abundance}
\keyword{OTU}
\keyword{otuTable}

