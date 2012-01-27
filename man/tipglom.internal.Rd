\name{tipglom.internal}
\alias{tipglom.internal}
\title{Internal function for tiplgom.}
\usage{
  tipglom.internal(tree, speciationMinLength)
}
\arguments{
  \item{tree}{An object of class \code{phylo}, or
  \code{phyloseq}}

  \item{speciationMinLength}{The minimum branch length that
  separates taxa. All tips of the tree separated by a
  cophenetic distance smaller than
  \code{speciationMinLength} will be agglomerated.}
}
\value{
  An object of class \code{phylo}, or \code{phyloseq}.
  Output class matches the class of \code{tree}.
}
\description{
  Internal function, users should use the S4 method
  \code{\link{tipglom}}. Tree can be a
  \code{\link{phyloseq-class}} that contains a phylogenetic
  tree, This is because \code{\link{merge_species}} can
  handle all the relevant objects, as can
  \code{\link{getTipDistMatrix}}. Create Non-trivial OTU
  table, by agglomerating nearby tips. tipglom.internal is
  called by the S4 \code{tipglom} methods. It is useful if
  a motivated user wants to see the internals of the
  implementation. By design it lacks explicit object
  handling. Use \code{\link{tipglom}} instead.
}
\seealso{
  tipglom
}
\keyword{internal}

