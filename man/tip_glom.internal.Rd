\name{tip_glom.internal}
\alias{tip_glom.internal}
\title{Internal function for tiplgom.}
\usage{
  tip_glom.internal(tree, speciationMinLength)
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
  \code{\link{tip_glom}}. Tree can be a
  \code{\link{phyloseq-class}} that contains a phylogenetic
  tree, This is because \code{\link{merge_taxa}} can handle
  all the relevant objects, as can
  \code{\link{getTipDistMatrix}}. Create Non-trivial OTU
  table, by agglomerating nearby tips. tip_glom.internal is
  called by the S4 \code{tip_glom} methods. It is useful if
  a motivated user wants to see the internals of the
  implementation. By design it lacks explicit object
  handling. Use \code{\link{tip_glom}} instead.
}
\seealso{
  tip_glom
}
\keyword{internal}

