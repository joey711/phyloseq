\name{sampleMap-class}
\alias{sampleMap-class}
\title{The S4 object in the phyloseq package that holds sample data as a data.frame.}
\description{
  Row indices represent samples, while column indices
  represent experimental categories, variables, etc. that
  describe the samples.
}
\details{
  \describe{ \item{nsamples}{ A single positive integer
  specifying the number of distinct samples.
  \code{nsamples} should not be assigned directly by user,
  but is assessed internally during initialization of
  \code{sampleMap} objects.}

  \item{.Data}{data-frame data, inherited from the
  data.frame class.}

  \item{row.names}{Also inherited from the data.frame
  class; it should contain the sample names.}

  \item{names}{Inherited from the data.frame class.}

  \item{nsamples}{A single positive integer, indicating the
  number of samples (rows) represented by an object.
  Consistent with otuTable objects.  should not be assigned
  directly by user, but is assessed internally during
  instantiation of \code{sampleMap} objects.} }
}

