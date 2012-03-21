\docType{methods}
\name{cca.phyloseq}
\alias{cca.phyloseq}
\alias{cca.phyloseq,formula-method}
\alias{cca.phyloseq,otuTable-method}
\alias{cca.phyloseq,phyloseq-method}
\alias{rda.phyloseq}
\alias{rda.phyloseq,formula-method}
\alias{rda.phyloseq,otuTable-method}
\alias{rda.phyloseq,phyloseq-method}
\title{Wrapper for \code{\link[vegan]{cca}} and \code{\link[vegan]{rda}}.}
\usage{
  cca.phyloseq(X, ...)

  rda.phyloseq(X, ...)
}
\arguments{
  \item{X}{(Required). A \code{\link{formula}}, specifying
  the input. No need to directly access components.
  \code{cca.phyloseq} understands where to find the
  abundance table and sample data. Alternatively, \code{X}
  can be an \code{\link{otuTable-class}} or
  \code{\link{phyloseq-class}} (without the \code{~}
  signifying a formula), in which case an unconstrained
  ordination is performed.}

  \item{...}{(Optional). E.g. \code{data=DF}, where
  \code{DF} is a \code{data.frame} containing information
  equivalent to a \code{sampleData} object / component.
  Only necessary if complex object does not already contain
  \code{sampleData} or you are keeping the data separate
  for some reason.}
}
\value{
  same output as \code{\link[vegan]{cca}} or
  \code{\link[vegan]{rda}}, respectively.
}
\description{
  A formula is main input to \code{\link[vegan]{cca}}. This
  complicates dispatch based on object signature. A new
  method with a separate name is defined instead.
}
\examples{
#
# data(GlobalPatterns)
# # For RDA, use thresholded-rank
# ex4  <- transformsamplecounts(GlobalPatterns, threshrankfun(500))
# # RDA
# modr <- rda.phyloseq(ex4 ~ SampleType)
# # CCA
# modc <- cca.phyloseq(GlobalPatterns ~ SampleType)
# plot_ordination(GlobalPatterns, modr, "biplot")
# plot_ordination(GlobalPatterns, modc, "biplot")
# # Perform unconstrained ordination
# mod1 <- cca.phyloseq(GlobalPatterns)
}
\seealso{
  \code{\link{plot_ordination}}, \code{\link[vegan]{rda}},
  \code{\link[vegan]{cca}}
}
\keyword{internal}

