\docType{methods}
\name{otuTable}
\alias{otutable}
\alias{otuTable}
\alias{otuTable,-method}
\alias{otuTable,data.frame-method}
\alias{otuTable,matrix-method}
\alias{otuTable,otuTable-method}
\alias{otuTable,phyloseq-method}
\title{Build or access otuTable objects.}
\usage{
  otuTable(object, speciesAreRows, errorIfNULL=TRUE)
}
\arguments{
  \item{object}{(Required). A phyloseq object.}

  \item{speciesAreRows}{(Conditionally optional). Ignored
  unless \code{object} is a matrix, in which case
  \code{speciesAreRows} is required.}

  \item{errorIfNULL}{(Optional). Logical. Should the
  accessor stop with an error if the slot is empty
  (\code{NULL})? Default \code{TRUE}. Ignored if
  \code{object} argument is a matrix (constructor invoked
  instead).}
}
\value{
  An otuTable object. It is either grabbed from the
  relevant slot if \code{object} is complex, or built anew
  if \code{object} is an integer matrix representing the
  species-abundance table. If \code{object} is a
  data.frame, then an attempt is made to coerce it to an
  integer matrix and instantiate an otuTable object.
}
\description{
  \code{otuTable()} is both a constructor and accessor
  method. When the first argument is a matrix, otuTable()
  will attempt to create and return an otuTable-class
  object, which further depends on whether or not
  speciesAreRows is provided as an additional argument.
  Alternatively, if the first argument is an object that
  contains an otuTable, including an otuTable-class object,
  then the corresponding otuTable is returned, as the
  component object by itself. This is a convenience wrapper
  on the more general \code{\link{access}} function
  specific for grabbing the otuTable of an object. It
  should work on both otuTable component objects, and
  higher-order classes that contain an otuTable slot.
}
\details{
  This is the main method suggested for constructing
  otuTable objects from an abundance matrix. It is also the
  suggested method for accessing subsetting the otuTable
  from a more complex object.
}
\examples{
#
# OTU1 <- otuTable(matrix(sample(0:5,250,TRUE),25,10), TRUE)
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

