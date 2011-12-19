\docType{methods}
\name{plot_phyloseq}
\alias{plot_phyloseq}
\alias{plot_phyloseq,phyloseq-method}
\title{Generic plot defaults for phyloseq.}
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
  non-empty slots. This is mainly for syntactic convenience
  and quick-plotting. See links below for some examples of
  available graphics tools available in the
  \code{\link{phyloseq-package}}.
}
\examples{
## data(ex1)
## plot_phyloseq(ex1)
}
\seealso{
  \code{\link{taxaplot}} \code{\link{plot_tree_phyloseq}}
  \code{\link{plot_ordination_phyloseq}}
  \code{\link{calcplot}} \code{\link{makenetwork}}
}

