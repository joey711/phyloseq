\name{tipglom.internal}
\alias{tipglom.internal}
\title{Internal function for tiplgom.}
\usage{
  tipglom.internal(tree, speciationMinLength)
}
\arguments{
  \item{tree}{An object of class \code{phylo},
  \code{phylo4}, \code{otuTree} or its superclasses.}

  \item{speciationMinLength}{The minimum branch length that
  separates taxa. All tips of the tree separated by a
  cophenetic distance smaller than
  \code{speciationMinLength} will be agglomerated.}
}
\value{
  An object of class \code{phylo}, \code{phylo4},
  \code{otuTree} or its superclasses. Output class matches
  the class of \code{tree}.
}
\description{
  Internal function, users should use the S4 method
  \code{\link{tipglom}}. Tree can be higher-order object
  that contains a phylogenetic tree, e.g. otuTree,
  otuTree4, etc. This is because \code{\link{mergespecies}}
  can handle all the relevant objects, as can
  \code{\link{getTipDistMatrix}}. Create Non-trivial OTU
  table, by agglomerating nearby tips. tipglom.internal is
  called by the S4 \code{tipglom} methods. It is useful if
  a motivated user wants to see the internals of the
  implementation. By design it lacks explicit object
  handling. Use \code{\link{tipglom}} instead.
}
\examples{
#
}
\seealso{
  tipglom
}

