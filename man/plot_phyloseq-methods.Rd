\docType{methods}
\name{plot_phyloseq}
\alias{plot_phyloseq}
\alias{plot_phyloseq,phyloseq-method}
\title{Generic plot defaults for phyloseq-objects.}
\usage{
  plot_phyloseq(physeq, ...)
}
\arguments{
  \item{physeq}{(Required). \code{\link{phyloseq-class}}.
  The actual plot type depends on the available (non-empty)
  component data types contained within.}

  \item{...}{(Optional). Additional parameters to be passed
  on to the respective specific plotting function. See
  below for different plotting functions that might be
  called by this generic plotting wrapper.}
}
\value{
  A plot is created. The nature and class of the plot
  depends on the \code{physeq} argument, specifically,
  which component data classes are present.
}
\description{
  The specific plot type is chosen according to available
  non-empty slots.
}
\examples{
## data(ex1)
## plot_phyloseq(ex1)
}

