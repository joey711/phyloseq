\name{access}
\alias{access}
\title{General slot accessor function for phyloseq-class.}
\usage{
  access(object, slot, errorIfNULL=FALSE)
}
\arguments{
  \item{object}{(Required). A phyloseq-class object.}

  \item{slot}{(Required). A character string indicating the
  slot (not data class) of the component data type that is
  desired.}

  \item{errorIfNULL}{(Optional). Logical. Should the
  accessor stop with an error if the slot is empty
  (\code{NULL})? Default \code{FALSE}.}
}
\value{
  Returns the component object specified by the argument
  \code{slot}. Returns NULL if slot does not exist. Returns
  \code{object} as-is if it is a component class that
  already matches the slot name.
}
\description{
  This function is used internally by many convenience
  accessors and in many functions/methods that need to
  access a particular type of component data. If something
  is wrong, or the slot is missing, the expected behavior
  is that this function will return NULL. Thus, the output
  can be tested by \code{\link{is.null}} as verification of
  the presence of a particular data component.
}
\examples{
#
## data(ex1)
## access(ex1, "taxTab")
## access(ex1, "tre")
## access(otuTable(ex1), "otuTable")
## # Should return NULL:
## access(otuTable(ex1), "sampleMap")
## access(otuTree(ex1), "sampleMap")
## access(otuSam(ex1), "tre")
}
\seealso{
  merge_phyloseq
}

