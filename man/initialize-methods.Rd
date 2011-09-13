\docType{methods}
\name{initialize}
\alias{initialize}
\alias{initialize,otuSam-method}
\alias{initialize,otuSamTax-method}
\alias{initialize,otuSamTaxTree-method}
\alias{initialize,otuSamTree-method}
\alias{initialize,otuTable-method}
\alias{initialize,otuTax-method}
\alias{initialize,otuTree-method}
\alias{initialize,sampleMap-method}
\alias{initialize,taxonomyTable-method}
\title{Create an otuTable object from an abundance matrix.}
\arguments{
  \item{.Object}{numeric (integer) matrix containing
  abundance data.}

  \item{speciesAreRows}{A single-element logical specifying
  the orientation of the abundance table.}

  \item{...}{Optional additional arguments. Ignored.}

  \item{.Object}{numeric (integer) matrix containing
  abundance data.}

  \item{.Object}{character matrix containing taxonomic
  classifiers, where each row is a different species / taxa
  in a high-throughput sequencing experiment.}

  \item{...}{Additional arguments. Ignored.}
}
\description{
  (first unnamed entry), and a logical specifying the
  orientation in a named entry for \code{speciesAreRows}.
  Other slots in the otuTable object are determined
  automatically from dimensions and dim-names. If no
  dim-names are provided, a default set is created with
  "sa" prefix for samples, and "sp" prefix for species /
  taxa.

  Create an sampleMap object from a required data.frame
  class. (first unnamed entry). Other slots in
  \code{\link{sampleMap}} object are determined
  automatically from dimensions and row.names. If no
  dim-names are provided, a default set is created with
  "sa" prefix for samples.

  If \code{.Object} has no dim-names, then taxonomic-level
  and species/taxa names are added as dummy default.
}
\details{
  This is an internal method. The suggested approach for
  creating an otuTable object from an abundance matrix is
  to use the \code{\link{otuTable}} method instead.
}
\examples{
## otuTable(matrix(sample(0:5,300,TRUE),30,10), speciesAreRows=FALSE)
## taxtab1 <- new("taxonomyTable", matrix("abc", 5, 5))
}
\seealso{
  \code{\link{otuTable}}

  \code{\link{sampleMap}}

  \code{\link{taxTab}}
}

