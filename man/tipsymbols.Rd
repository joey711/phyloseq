\name{tipsymbols}
\alias{tipsymbols}
\title{Annotate tips on a tree.}
\usage{
  tipsymbols(tip, adj = c(0.5, 0.5), ...)
}
\arguments{
  \item{tip}{An integer specifying the tip ID in a tree
  that for which the base plot has already been generated
  and is still available to \code{R}.}

  \item{adj}{A 2 element numeric vector specifying a
  position adjustment.}

  \item{...}{Additional plotting parameters that are passed
  to \code{\link{points}} or \code{\link{text}} in the R
  base graphics.}
}
\value{
  No objects returned. Symbol or text is plotted on the
  available graphic device.
}
\description{
  There were some unexpected behavior from the tiplabels
  function in ape. These functions are intended to act as
  simplified versions that act as a convenience wrapper for
  \code{points()} or \code{text()} functions, respectively,
  but where the tip coordinates are specified by giving the
  tip ID (integer) as input.
}
\examples{
#
## tipsymbols(1, pch=19)
## tiptext(1, labels="my.label")
}
\seealso{
  tiplabels points text
}

