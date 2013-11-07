\name{metaMDS}
\alias{metaMDS}
\title{S3 class placeholder definition (list) for metaMDS}
\description{
  The ape package does export a version of its
  \code{\link[vegan]{metaMDS}}-class, partly because it is
  not really defined formally anywhere. Instead, it is an
  S3 class extended from the base class, \code{\link{list}}
  -- this is a very common and easy approach -- and proper
  behavior of any method taking an instance of this class
  requires exact naming conventions for element names of
  the list components. The phyloseq package does not
  provide any validity checks that a given phylo instance
  is valid (conforms to the conventions in the ape
  package)... yet. If problems arise, this might be
  considered, and they could be defined judiciously and
  within phyloseq.
}
\seealso{
  \code{\link[vegan]{metaMDS}}
}
\keyword{internal}

