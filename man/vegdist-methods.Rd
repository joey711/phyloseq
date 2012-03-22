\docType{methods}
\name{vegdist}
\alias{vegdist}
\alias{vegdist,otuTable-method}
\alias{vegdist,phyloseq-method}
\title{\code{\link[vegan]{vegdist}} wrapper for phyloseq classes}
\description{
  Trivially-extended S4 method from the
  \code{\link[vegan]{vegdist}} function, such that S4
  classes from the \code{\link{phyloseq-package}} are
  properly handled / accessed. All parameters passed on to
  \code{\link[vegan]{vegdist}} verbatim.
}
\examples{
# data(esophagus)
# vegdist(esophagus, "jaccard")
}
\seealso{
  \code{\link[vegan]{vegdist}}
}
\keyword{internal}

