\docType{methods}
\name{otuTable}
\alias{otutable}
\alias{otuTable}
\alias{otuTable,-method}
\alias{otuTable,data.frame-method}
\alias{otuTable,matrix-method}
\alias{otuTable,otuTable-method}
\alias{otuTable,phyloseqFather-method}
\title{Build or access otuTable objects.}
\arguments{
  \item{object}{An object among the set of classes defined
  by the phyloseq package that contain an otuTable.}

  \item{...}{(optional) ignored unless \code{object} is a
  matrix, in which case \code{speciesAreRows} is required
  as a named argument.}

  \item{object}{An otuTable-class object. The abundance
  table.}

  \item{object}{Any \emph{phyloseq}-package object that
  contains an otuTable.}

  \item{object}{An abundance table in the form of an
  integer matrix. # # \code{speciesAreRows} must be
  specified in \code{...}-argument . # #}

  \item{...}{The additional named argument
  \code{speciesAreRows} must be # # provided. A logical.}
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

  # # Return the otuTable object as-is. # # # #

  # # Return the otuTable component from a H.O. object. # #
  # #

  # # Instantiate an otuTable from a raw abundance matrix.
  # # # #

  # # Convert to matrix, then instantiate otuTable. # #
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

