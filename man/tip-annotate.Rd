\name{tipsymbols}
\alias{tipsymbols}
\alias{tiptext}
\title{Annotate tips on a tree with symbols or text.}
\usage{
  tipsymbols(tip, adj=c(0.5, 0.5), ...)

  tiptext(tip, adj = c(0.5, 0.5), ...)
}
\arguments{
  \item{tip}{An integer specifying the tip ID in a tree
  that for which the base plot has already been generated
  and is still available to \code{R}.}

  \item{adj}{A 2 element numeric vector specifying a
  position adjustment.}

  \item{...}{Additional plotting parameters that are passed
  to \code{\link{points}} or \code{\link{text}} in the R
  base graphics. Again, for \code{tiptext()}, make sure to
  include a \code{labels=} argument.}
}
\value{
  No objects returned. Symbol or text is plotted on the
  available graphic device.
}
\description{
  There were some unexpected behavior from the
  \code{\link[ape]{tiplabels}} function in ape. These
  functions are intended to act as simplified versions that
  act as a convenience wrapper for \code{points()} or
  \code{text()} functions, respectively, but where the tip
  coordinates are specified by giving the tip ID (integer)
  as input. For \code{tiptext()}, make sure to include a
  \code{labels=} argument, which will be passed on to
  \code{\link[graphics]{text}}.
}
\examples{
#
## data(GlobalPatterns)
## # for reproducibility
## set.seed(711)
## ex2 <- prune_species(sample(species.names(GlobalPatterns), 50), GlobalPatterns)
## plot( tre(ex2) )
## tipsymbols(pch=19)
## tipsymbols(1, pch=22, cex=3, col="red", bg="blue")
## tiptext(2, labels="my.label")
}
\seealso{
  \code{\link[ape]{tiplabels}},
  \code{\link[graphics]{points}},
  \code{\link[graphics]{text}}
}

