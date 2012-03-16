\docType{methods}
\name{sample.names}
\alias{sample.names}
\alias{sample.names,ANY-method}
\alias{sample.names,otuTable-method}
\alias{sample.names,phyloseq-method}
\alias{sample.names,sampleData-method}
\alias{sampleNames}
\title{Get sample names.}
\usage{
  sample.names(physeq)
}
\arguments{
  \item{physeq}{(Required). A \code{\link{phyloseq-class}},
  \code{\link{sampleData}}, or
  \code{\link{otuTable-class}}.}
}
\value{
  A character vector. The names of the samples in
  \code{physeq}.
}
\description{
  Get sample names.
}
\examples{
#
# # From "picante" package
# data(GlobalPatterns)
# sample.names(GlobalPatterns)
}
\seealso{
  \code{\link{species.names}}, \code{\link{nsamples}}
}

