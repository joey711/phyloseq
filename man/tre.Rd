\name{tre}
\alias{tre}
\title{Access phylo-class objects.}
\arguments{
  \item{object}{An object among the set of classes defined
  by the phyloseq package that contain a phylogenetic
  tree.}
}
\value{
  A phylo object. It is grabbed from the tre slot
}
\description{
  \code{tre()} is an accessor method. This is the main
  method suggested for accessing/subsetting the
  phylogenetic tree (phylo class) from a more complex
  object. \code{tre} is the slot name that holds the
  phylo-class object in a multi-component phyloseq-package
  object.
}
\details{
  Note that the tip.labels should be named to match the
  \code{species.names} of the other objects to which it is
  paired. The initialization methods for more complex
  objects automatically check for exact agreement in the
  set of species described by the phlyogenetic tree and the
  other components (taxonomyTable, otuTable). They also
  trim accordingly. Thus, the tip.labels in a phylo object
  must be named to match the \code{species.names} of the
  other objects to which it will ultimately be paired.
}
\seealso{
  otuTable sampleMap taxTab phyloseq
}
\keyword{OTU}
\keyword{phylo}
\keyword{tre}
\keyword{tree}

