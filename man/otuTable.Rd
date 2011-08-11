\name{otuTable}
\alias{otutable}
\alias{otuTable}
\title{The S4 class that holds taxa-abundance information.}
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
  Because orientation of these tables can vary by method,
  the orientation is defined explicitly in the
  \code{speciesAreRows} slot (a logical). The
  \code{otuTable} class inherits the \code{\link{matrix}}
  class to store abundance values. Various standard subset
  and assignment nomenclature has been extended to apply to
  the \code{otuTable} class, including square-bracket,
  \code{\link{t}}, etc.

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
  at-slot: speciesAreRows A single logical specifying the
  orientation of the abundance table.

  at-slot: nspecies A single positive integer specifying
  the number of distinct species/taxa.

  at-slot: nsamples A single positive integer specifying
  the number of distinct samples. Both \code{nsamples} and
  \code{nspecies} should not be assigned directly by user,
  but are assessed internally during initialization of
  \code{otuTable} objects.

  at-slot: species.names A character vector of the names of
  the taxa / species, so that they can be accessed,
  modified, subsetted, without the need for the user to
  consider the orientation of the abundance table.

  at-slot: sample.names A character vector of the names of
  the samples, so that they can be accessed, modified,
  subsetted, without the need for the user to consider the
  orientation of the abundance table. Both
  \code{sample.names} and \code{species.names} can be
  assigned, e.g. \code{species.names(object)<-}, however,
  this has the effect of internally replacing the otuTable
  to be a subset that matches the new name vector, rather
  than replacing the index names in-place.

  at-slot: .Data This slot is inherited from the
  \code{\link{matrix}} class.

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
\keyword{internal}
\keyword{OTU}
\keyword{otuTable}

