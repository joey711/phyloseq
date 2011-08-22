\name{tipglom}
\alias{tipglom}
\title{Agglomerate nearby tips in tree.}
\arguments{
  \item{tree}{An object of class \code{otuTree} or its
  superclasses, in which case the OTU argument can be
  ommitted. If, alternatively, \code{tree} is a
  \code{phylo} or \code{phylo4} object, then \code{OTU} is
  required.}

  \item{OTU}{optional. Ignored if \code{tree} is an
  \code{otuTree} object or one of its superclasses.
  Required if \code{tree} is a \code{phylo} or
  \code{phylo4} object. \code{OTU} must be an
  \code{otuTable} object or a 2-column table that indicates
  the sample from which each sequence originated.}

  \item{speciationMinLength}{The minimum branch length that
  separates taxa. All tips of the tree separated by a
  cophenetic distance smaller than
  \code{speciationMinLength} will be agglomerated. Default
  is 0.02}
}
\value{
  An object of class \code{otuTree} or its superclasses.
  Output class matches the class of \code{tree}, unless it
  is a \code{phylo}/\code{phylo4} object, in which case
  \code{tipglom} returns an \code{otuTree} object.
}
\description{
  tipglom determines those groups of tips in the tree that
  are closer than the threshold provided as
  \code{speciationMinLength}. That is, all tips of the tree
  separated by a cophenetic distance smaller than
  \code{speciationMinLength} will be agglomerated using
  \code{mergespecies}. Can be used to create non-trivial
  OTU table, by agglomerating nearby tips.
}
\examples{
#
}

