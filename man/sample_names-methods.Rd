\docType{methods}
\name{sample_names}
\alias{sample_names}
\alias{sample_names,ANY-method}
\alias{sample_names,otu_table-method}
\alias{sample_names,phyloseq-method}
\alias{sample_names,sample_data-method}
\alias{sample.names}
\alias{sampleNames}
\title{Get sample names.}
\usage{
  sample_names(physeq)

  sampleNames(physeq)

  sample.names(physeq)
}
\arguments{
  \item{physeq}{(Required). A \code{\link{phyloseq-class}},
  \code{\link{sample_data}}, or
  \code{\link{otu_table-class}}.}
}
\value{
  A character vector. The names of the samples in
  \code{physeq}.
}
\description{
  Get sample names.
}
\examples{
# data(GlobalPatterns)
# sample_names(GlobalPatterns)
}
\seealso{
  \code{\link{taxa_names}}, \code{\link{nsamples}}
}

