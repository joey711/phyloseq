\docType{methods}
\name{sampleMap}
\alias{samplemap}
\alias{sampleMap}
\alias{sampleMap,ANY-method}
\alias{sampleMap,data.frame-method}
\title{Build or access sampleMap objects.}
\arguments{
  \item{object}{An object among the set of classes defined
  by the phyloseq package that contain sampleMap.}
}
\value{
  A sampleMap object. It is either grabbed from the
  relevant slot if \code{object} is complex, or built anew
  if \code{object} is a data.frame representing the sample
  covariates of an experiment.
}
\description{
  \code{sampleMap()} is both a constructor and accessor
  method. When the argument is a data.frame, sampleMap()
  will attempt to create and return an sampleMap-class
  object. In this case, the rows should be named to match
  the \code{sample.names} of the other objects to which it
  will ultimately be paired. Alternatively, if the argument
  is an object that contains a sampleMap, including a
  sampleMap-class object, then the corresponding sampleMap
  is returned, as-is.
}
\details{
  This is the main method suggested for constructing
  sampleMap objects from a data.frame of sample covariates.
  Each row represents a different sample. It is also the
  suggested method for accessing/subsetting the sampleMap
  from a more complex object.
}
\examples{
#
# OTU1 <- otuTable(matrix(sample(0:5,250,TRUE),25,10), speciesAreRows=TRUE)
# tax1 <- taxTab(matrix("abc", 30, 8))
# map1 <- data.frame( matrix(sample(0:3,250,TRUE),25,10),
#   matrix(sample(c("a","b","c"),150,TRUE), 25, 6) )
# map1 <- sampleMap(map1)
# ex1 <- phyloseq(OTU1, map1, tax1)
# sampleMap(ex1)
}
\seealso{
  otuTable taxTab tre phyloseq sampleMap-class
}

